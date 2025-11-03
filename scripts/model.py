# scripts/model.py
# --- NEW VERSION ---
# This model file replaces the simple GNN with the ENINet (EquiThreeBody) architecture.
#
# IMPORTANT: This code requires the 'eninet' package and 'dgl' to be installed.
# Ensure you have run 'pip install -e .' from the ENINet repo and 'pip install dgl'.

from __future__ import annotations

import torch
import torch.nn as nn
import dgl
import mdtraj as md
import numpy as np

# --- Imports from ENINet package ---
# These are required for the EquiThreeBody model.
try:
    from eninet.layer import (
        CosineCutoff,
        GaussianRBF,
        ThreeBodyEquiGraphConvSimple,
        TwoBodyEquiGraphConv
    )
except ImportError:
    print("ERROR: Could not import from 'eninet.layer'.")
    print("Please make sure you have installed the ENINet package by running 'pip install -e .' inside its repository.")
    raise

from dgl import DGLGraph

# ---
# 1. The Core ENINet (EquiThreeBody) Engine
#    (This is the class from the ENINet repo's model.py)
# ---

class EquiThreeBody(torch.nn.Module):
    """
    The EquiThreeBody (ENINet) GNN engine.
    This module performs the main equivariant message passing.
    """
    def __init__(
        self,
        n_elements: int,
        g_feat_dim: int = 128,
        lg_feat_dim: int = 16,
        n_interactions: int = 3,
        n_rbf: int = 20,
        cutoff=5.0,  # Note: This cutoff is in Angstroms
        activation: str = "silu",
        g_aggregation: str = "sum",
        lg_aggreation: str = "sum",
        use_linegraph: bool = True,
    ) -> None:
        super().__init__()

        self.g_feat_dim = g_feat_dim
        self.lg_feat_dim = lg_feat_dim
        self.n_interactions = n_interactions
        self.cutoff = cutoff
        self.use_linegraph = use_linegraph

        # Cutoff and radial basis functions (RBF)
        self.twobody_cutoff_fn = CosineCutoff(cutoff)
        self.threebody_cutoff_fn = CosineCutoff(cutoff * 2)
        self.twobody_dist_rb = GaussianRBF(0, cutoff, n_rbf)
        self.threebody_dist_rb = GaussianRBF(0, cutoff * 2, n_rbf)

        # Embeddings for atoms, bonds, and triplets
        self.n_embedding = nn.Embedding(n_elements, g_feat_dim)
        self.e_embedding = nn.Linear(n_rbf, g_feat_dim)
        self.t_embedding = nn.Linear(n_rbf, lg_feat_dim)

        # Interaction blocks
        self.twobody_conv = nn.ModuleList(
            [
                TwoBodyEquiGraphConv(
                    feat_dim=g_feat_dim,
                    n_rbf=n_rbf,
                    activation=activation,
                    aggregation=g_aggregation,
                    cutoff=cutoff,
                )
                for _ in range(n_interactions)
            ]
        )

        if use_linegraph:
            self.threebody_conv = nn.ModuleList(
                [
                    ThreeBodyEquiGraphConvSimple(
                        g_feat_dim=g_feat_dim,
                        lg_feat_dim=lg_feat_dim,
                        n_rbf=n_rbf,
                        activation=activation,
                        aggregation=lg_aggreation,
                        cutoff=cutoff * 2,
                    )
                    for _ in range(n_interactions)
                ]
            )

        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_uniform_(self.n_embedding.weight)
        nn.init.xavier_uniform_(self.e_embedding.weight)
        nn.init.zeros_(self.e_embedding.bias)
        nn.init.xavier_uniform_(self.t_embedding.weight)
        nn.init.zeros_(self.t_embedding.bias)

    def _compute_neighbor_features(self, g):
        pos_j = g.ndata["pos"][g.edges()[1]]
        if "pbc_offshift" in g.edata:
            pos_j += g.edata["pbc_offshift"]
        pos_i = g.ndata["pos"][g.edges()[0]]
        vctr_ij = pos_j - pos_i
        dist_ij = vctr_ij.norm(dim=1, keepdim=True)

        # Avoid division by zero for zero-distance pairs (if any)
        dist_ij = dist_ij.clamp(min=1e-6)

        g.edata["vctr_norm"] = (vctr_ij / dist_ij) * self.cutoff
        g.edata["dist"] = dist_ij
        return g

    def _compute_triplet_features(self, g, l_g):
        pos_k = g.ndata["pos"][g.find_edges(l_g.edges()[1])[1]]
        if "pbc_offshift" in g.edata:
            pos_k += (
                g.edata["pbc_offshift"][l_g.edges()[0]]
                + g.edata["pbc_offshift"][l_g.edges()[1]]
            )
        pos_j = g.ndata["pos"][g.find_edges(l_g.edges()[0])[0]]
        vctr_jk = pos_k - pos_j
        dist_jk = vctr_jk.norm(dim=1, keepdim=True)

        # Avoid division by zero
        dist_jk = dist_jk.clamp(min=1e-6)

        l_g.edata["vctr_norm"] = vctr_jk / dist_jk
        l_g.edata["dist"] = dist_jk
        return l_g

    def _init_atoms(self, g):
        node_s = self.n_embedding(g.ndata["node_type"])[:, None]
        node_v = torch.zeros((node_s.size(0), 3, node_s.size(2)), device=node_s.device)
        return node_s, node_v

    def _init_bonds(self, g):
        dist_ij = g.edata["dist"]
        rb_ij = self.twobody_dist_rb(dist_ij)
        fcut_ij = self.twobody_cutoff_fn(dist_ij)

        edge_s = self.e_embedding(rb_ij) * fcut_ij[..., None]
        edge_v = (
            g.edata["vctr_norm"][..., None].expand(-1, -1, self.g_feat_dim)
            * fcut_ij[..., None]
        )
        return edge_s, edge_v

    def _init_triplets(self, l_g):
        dist_jk = l_g.edata["dist"]
        rb_jk = self.threebody_dist_rb(dist_jk)
        fcut_jk = self.threebody_cutoff_fn(dist_jk)

        triplet_s = self.t_embedding(rb_jk) * fcut_jk[..., None]
        triplet_v = (
            l_g.edata["vctr_norm"][..., None].expand(-1, -1, self.lg_feat_dim)
            * fcut_jk[..., None]
        )
        return triplet_s, triplet_v

    def forward(self, g: DGLGraph, l_g: DGLGraph) -> DGLGraph:
        g = self._compute_neighbor_features(g)
        node_s, node_v = self._init_atoms(g)
        edge_s, edge_v = self._init_bonds(g)

        if self.use_linegraph:
            l_g = self._compute_triplet_features(g, l_g)
            triplet_s, triplet_v = self._init_triplets(l_g)

        for i in range(self.n_interactions):
            if self.use_linegraph:
                edge_s, edge_v, triplet_s, triplet_v = self.threebody_conv[i](
                    l_g, edge_s, edge_v, triplet_s, triplet_v
                )

            node_s, node_v, edge_s, edge_v = self.twobody_conv[i](
                g, node_s, node_v, edge_s, edge_v
            )

        # Squeeze scalar features for readout
        g.ndata["s"] = node_s.squeeze(1)
        # Permute vector features for readout
        g.ndata["v"] = node_v.permute(0, 2, 1) # Shape: [n_atoms, g_feat_dim, 3]
        return g


# ---
# 2. The CommittorNet Wrapper
#    (This is the class your script will import and use)
# ---

class CommittorNet(nn.Module):
    """
    A wrapper for EquiThreeBody to predict the committor probability.

    This model takes a DGLGraph (g) and its line graph (l_g),
    runs the equivariant GNN, and then uses a readout head to output a
    single value between 0 and 1.
    """
    def __init__(self,
                 num_atom_types=100,
                 embedding_dim=64,
                 n_layers=3,
                 lg_embedding_dim=16,
                 cutoff_A=5.0): # Cutoff is now in Angstroms
        super(CommittorNet, self).__init__()

        # --- GNN Engine ---
        self.equi_three_body = EquiThreeBody(
            n_elements=num_atom_types,
            g_feat_dim=embedding_dim,
            lg_feat_dim=lg_embedding_dim,
            n_interactions=n_layers,
            n_rbf=20, # Default from EquiThreeBody
            cutoff=cutoff_A, # Cutoff in Angstroms
            use_linegraph=True
        )

        # --- Readout Head ---
        # This part takes the final atom features (scalar and vector)
        # and pools them into a single, invariant graph-level prediction.

        # Input to readout will be:
        # 1. Pooled scalar features (dim: embedding_dim)
        # 2. Norm of pooled vector features (dim: embedding_dim)
        readout_input_dim = embedding_dim + embedding_dim

        self.readout = nn.Sequential(
            nn.Linear(readout_input_dim, 32),
            nn.SiLU(),  # SiLU is a common activation in these models
            nn.Linear(32, 1),
            nn.Sigmoid()
        )

    def forward(self, g, l_g):
        # 1. Run the equivariant GNN engine
        # This adds the final features 's' and 'v' to the graph 'g'
        g = self.equi_three_body(g, l_g)

        # 2. Pool the node features to get graph-level features
        # We use 'mean' pooling, which is invariant

        # Pool scalar features (invariant)
        # g.ndata['s'] has shape [n_atoms, embedding_dim]
        pooled_s = dgl.readout_nodes(g, 's', op='mean')

        # Pool vector features (equivariant)
        # g.ndata['v'] has shape [n_atoms, embedding_dim, 3]
        pooled_v = dgl.readout_nodes(g, 'v', op='mean')

        # 3. Convert pooled vector features to an invariant representation
        # We do this by taking the norm of the vectors for each feature channel.
        # This results in a tensor of shape [batch_size, embedding_dim]
        pooled_v_norm = torch.norm(pooled_v, p=2, dim=-1) # Norm across the 3 (x,y,z) dimensions

        # 4. Concatenate invariant features and make prediction
        # Shape: [batch_size, embedding_dim + embedding_dim]
        final_graph_features = torch.cat([pooled_s, pooled_v_norm], dim=1)

        # Return the final committor probability
        return self.readout(final_graph_features)


# ---
# 3. The New Data Preprocessing Function
#    (This function is also imported by your script)
# ---

# Define the cutoff in Angstroms, matching the model's default
GRAPH_CUTOFF_A = 5.0

def frame_to_torch_graph(frame: md.Trajectory, device: torch.device):
    """
    Converts a single mdtraj frame to a DGLGraph (g) and its line graph (l_g).
    
    CRITICAL FIX: Selects only protein atoms to avoid OOM errors from massive
    solvated systems (water/ions can be ~490k atoms).

    Args:
        frame: An mdtraj.Trajectory object with a single frame.
        device: The torch device to send the tensors to.

    Returns:
        A tuple of (g, l_g) for the CommittorNet.
    """
    # --- FIX: Select only protein atoms (exclude water, ions, lipids) ---
    # This reduces the graph from ~493k atoms to ~5k-10k protein atoms
    protein_indices = frame.topology.select('protein')
    
    if len(protein_indices) == 0:
        raise ValueError("No protein atoms found in trajectory frame!")
    
    # Create a sub-trajectory with only protein atoms
    protein_frame = frame.atom_slice(protein_indices)
    
    # 1. Get atom types and coordinates for PROTEIN ONLY
    atom_types = torch.tensor(
        [atom.element.atomic_number for atom in protein_frame.topology.atoms],
        dtype=torch.long,
        device=device,
    )

    # --- CRITICAL ---
    # Convert coordinates from nm (mdtraj default) to Angstroms (model default)
    coords_A = torch.tensor(protein_frame.xyz[0] * 10.0, dtype=torch.float32, device=device)
    num_atoms = len(atom_types)
    
    print(f"[Graph Construction] Using {num_atoms} protein atoms (reduced from {frame.n_atoms} total atoms)")

    # 2. Create the molecular graph (g) using KDTree for efficiency
    # For protein-only systems (~5k-10k atoms), we use scipy KDTree to avoid
    # the O(N^2) memory cost of torch.cdist in dgl.radius_graph
    from scipy.spatial import cKDTree
    
    # Build KDTree on CPU for neighbor search
    coords_np = coords_A.cpu().numpy()
    tree = cKDTree(coords_np)
    
    # Find all pairs within cutoff distance
    pairs = tree.query_pairs(r=GRAPH_CUTOFF_A, output_type='ndarray')
    
    if len(pairs) == 0:
        # Fallback: create a minimal graph with self-loops if no edges found
        print(f"[WARNING] No edges found within {GRAPH_CUTOFF_A}Å cutoff. Creating minimal graph.")
        src = torch.arange(num_atoms, device=device)
        dst = torch.arange(num_atoms, device=device)
    else:
        # Create bidirectional edges (i->j and j->i)
        src = np.concatenate([pairs[:, 0], pairs[:, 1]])
        dst = np.concatenate([pairs[:, 1], pairs[:, 0]])
        src = torch.from_numpy(src).long().to(device)
        dst = torch.from_numpy(dst).long().to(device)
    
    # Create DGL graph from edge list
    g = dgl.graph((src, dst), num_nodes=num_atoms).to(device)
    g.ndata['node_type'] = atom_types
    g.ndata['pos'] = coords_A

    # 3. Create the line graph (l_g)
    # This graph represents 3-body interactions (angles)
    # Its nodes are the edges of 'g', and its edges are connections
    # between edges in 'g' that share a common atom.
    l_g = dgl.line_graph(g, backtracking=False)

    print(f"[Graph Construction] Created graph: {g.num_nodes()} nodes, {g.num_edges()} edges")
    print(f"[Graph Construction] Line graph: {l_g.num_nodes()} nodes, {l_g.num_edges()} edges")

    # Return the two graphs, ready to be passed to the model
    return g, l_g
