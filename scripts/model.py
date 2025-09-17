# scripts/model.py
import torch
import torch.nn as nn


class EGNNLayer(nn.Module):
    """A single layer of a simple Equivariant Graph Neural Network."""

    def __init__(self, in_features, hidden_features, out_features):
        super().__init__()
        self.message_net = nn.Sequential(
            nn.Linear(in_features * 2 + 1, hidden_features),
            nn.SiLU(),
            nn.Linear(hidden_features, hidden_features),
            nn.SiLU(),
        )
        self.update_net = nn.Sequential(
            nn.Linear(in_features + hidden_features, out_features),
            nn.SiLU(),
            nn.Linear(out_features, out_features),
        )

    def forward(self, h, coords, edge_index):
        # h: node features (scalar)
        # coords: node coordinates (vector)
        # edge_index: graph connectivity

        sender, receiver = edge_index
        rel_coords = coords[sender] - coords[receiver]
        dist = torch.norm(rel_coords, dim=1, keepdim=True)

        # Create messages
        message_input = torch.cat([h[sender], h[receiver], dist], dim=1)
        messages = self.message_net(message_input)

        # Aggregate messages
        agg_messages = torch.zeros_like(h)
        agg_messages.index_add_(0, receiver, messages)

        # Update node features
        update_input = torch.cat([h, agg_messages], dim=1)
        h_new = h + self.update_net(update_input)

        return h_new


class CommittorNet(nn.Module):
    """The main GNN to predict the committor."""

    def __init__(self, in_features=128, hidden_features=128, num_layers=3):
        super().__init__()
        self.embedding = nn.Embedding(100, in_features)  # For atom types

        self.layers = nn.ModuleList(
            [EGNNLayer(in_features, hidden_features, in_features) for _ in range(num_layers)]
        )

        self.output_head = nn.Sequential(
            nn.Linear(in_features, hidden_features), nn.SiLU(), nn.Linear(hidden_features, 1)
        )

    def forward(self, atom_types, coords, edge_index):
        h = self.embedding(atom_types)

        for layer in self.layers:
            h = layer(h, coords, edge_index)

        # The output is a single logit value per graph
        graph_embedding = h.mean(dim=0)  # Average all node features
        logit = self.output_head(graph_embedding)

        # Apply sigmoid to get a probability between 0 and 1
        return torch.sigmoid(logit)
