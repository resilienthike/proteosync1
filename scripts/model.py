# scripts/model.py
import torch
import torch.nn as nn
from scipy.spatial import KDTree
import numpy as np
import mdtraj as md

class CommittorNet(nn.Module):
    """
    A Graph Isomorphism Network for predicting the committor probability.
    This model takes atomic graph data and outputs a single value between 0 and 1.
    """
    def __init__(self, num_atom_types=100, embedding_dim=64, n_layers=3):
        super(CommittorNet, self).__init__()
        self.embedding = nn.Embedding(num_atom_types, embedding_dim)
        
        # GNN layers for message passing
        self.gnn_layers = nn.ModuleList([
            nn.Sequential(
                nn.Linear(embedding_dim, embedding_dim),
                nn.ReLU(),
                nn.Linear(embedding_dim, embedding_dim)
            ) for _ in range(n_layers)
        ])
        
        # Readout network to produce a final prediction
        self.readout = nn.Sequential(
            nn.Linear(embedding_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )

    def forward(self, atom_types, coords, edge_index):
        # Note: This is a simplified GNN. For a full implementation inspired by the
        # ENINet paper, you would incorporate the coordinates (coords) in an
        # equivariant way. For now, we focus on the core AIMMD logic.
        h = self.embedding(atom_types)
        
        for layer in self.gnn_layers:
            # Aggregate information from neighbors
            row, col = edge_index
            agg = torch.zeros_like(h)
            
            # Sum messages from neighbors
            agg = agg.index_add(0, row, h[col])
            
            # Normalize by degree to get the mean
            degree = torch.zeros(h.size(0), dtype=torch.float32, device=h.device)
            degree = degree.index_add_(0, row, torch.ones_like(row, dtype=torch.float32))
            agg = agg / degree.unsqueeze(-1).clamp(min=1)

            # Update node embeddings (residual connection)
            h = h + layer(agg)

        # Pool node features to get a single graph-level representation
        graph_embedding = h.mean(dim=0)
        return self.readout(graph_embedding)

def frame_to_torch_graph(frame: md.Trajectory, device: torch.device):
    """
    Converts a single mdtraj frame to a graph data format for our GNN.
    
    Args:
        frame: An mdtraj.Trajectory object with a single frame.
        device: The torch device to send the tensors to.

    Returns:
        A tuple of (atom_types, coordinates, edge_index).
    """
    atom_types = torch.tensor([atom.element.atomic_number for atom in frame.topology.atoms], dtype=torch.long, device=device)
    coords = torch.tensor(frame.xyz[0], dtype=torch.float32, device=device)
    
    # Build edges based on a distance cutoff
    cutoff_nm = 0.5  # 5 Angstroms
    kdtree = KDTree(coords.cpu().numpy())
    edge_pairs = kdtree.query_pairs(r=cutoff_nm)
    
    if not edge_pairs:
        return atom_types, coords, torch.empty((2, 0), dtype=torch.long, device=device)

    # GNNs require an undirected graph, so we add edges in both directions
    adj = np.array(list(edge_pairs), dtype=np.int64)
    adj_rev = np.fliplr(adj)
    edge_index = torch.from_numpy(np.vstack((adj, adj_rev))).t().contiguous().to(device)
    
    return atom_types, coords, edge_index
