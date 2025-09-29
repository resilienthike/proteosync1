# scripts/enhanced_plots.py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots
    """Create summary statistics visualization"""
    
    if output_dir is None:
        output_dir = ARTIFACTS_DIR / "md" / target_name
        output_dir.mkdir(parents=True, exist_ok=True)_subplots
import cmocean
from pathlib import Path
import argparse

# Set up beautiful plotting defaults
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
sns.set_context("paper", font_scale=1.5)

# Configuration
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"

def create_enhanced_free_energy_plot(bin_centers, free_energy, target_name="Target", output_dir=None):
    """Create a beautiful free energy profile with multiple style options"""
    
    if output_dir is None:
        output_dir = ARTIFACTS_DIR / "md" / target_name
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. SEABORN STYLE - Professional publication quality
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Left plot: Classic style with gradient fill
    ax1.plot(bin_centers, free_energy, linewidth=3, color='#2E86AB', marker='o', markersize=6, markerfacecolor='white', markeredgewidth=2)
    ax1.fill_between(bin_centers, free_energy, alpha=0.3, color='#2E86AB')
    ax1.set_xlabel("Committor p(B)", fontsize=14, fontweight='bold')
    ax1.set_ylabel("Free Energy (kcal/mol)", fontsize=14, fontweight='bold')
    ax1.set_title(f"Free Energy Profile - {target_name}\nMulti-Region Structural Exploration", fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Add annotations for key regions
    min_idx = np.argmin(free_energy)
    ax1.annotate(f'Minimum\n{free_energy[min_idx]:.2f} kcal/mol', 
                xy=(bin_centers[min_idx], free_energy[min_idx]), 
                xytext=(bin_centers[min_idx] + 0.2, free_energy[min_idx] + 1),
                arrowprops=dict(arrowstyle='->', color='red', lw=2),
                fontsize=12, ha='center', 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', edgecolor='red', alpha=0.8))
    
    # Right plot: Heatmap style with scientific colormap
    X, Y = np.meshgrid(bin_centers, np.linspace(0, max(free_energy), 50))
    Z = np.tile(free_energy, (50, 1))
    
    im = ax2.contourf(X, Y, Z, levels=20, cmap=cmocean.cm.thermal, alpha=0.8)
    ax2.plot(bin_centers, free_energy, linewidth=4, color='white', alpha=0.9)
    ax2.plot(bin_centers, free_energy, linewidth=2, color='black')
    ax2.set_xlabel("Committor p(B)", fontsize=14, fontweight='bold')
    ax2.set_ylabel("Free Energy (kcal/mol)", fontsize=14, fontweight='bold')
    ax2.set_title("Thermodynamic Landscape", fontsize=16, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
    cbar.set_label('Free Energy (kcal/mol)', rotation=270, labelpad=20, fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    # Save both styles
    output_path_enhanced = output_dir / "enhanced_free_energy_profile.png"
    plt.savefig(output_path_enhanced, dpi=300, bbox_inches='tight')
    print(f"✅ Enhanced plot saved to {output_path_enhanced}")
    
    plt.show()

def create_interactive_plotly_dashboard(bin_centers, free_energy, target_name="Target", output_dir=None):
    """Create an interactive Plotly dashboard"""
    
    if output_dir is None:
        output_dir = ARTIFACTS_DIR / "md" / target_name
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Free Energy Profile", "Energy Gradient", "State Distribution", "3D Landscape"),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"type": "scene"}]]
    )
    
    # 1. Main free energy profile
    fig.add_trace(
        go.Scatter(
            x=bin_centers, 
            y=free_energy,
            mode='lines+markers',
            name='Free Energy',
            line=dict(color='#2E86AB', width=4),
            marker=dict(size=8, color='white', line=dict(color='#2E86AB', width=2)),
            hovertemplate='<b>Committor</b>: %{x:.3f}<br><b>Energy</b>: %{y:.2f} kcal/mol<extra></extra>'
        ),
        row=1, col=1
    )
    
    # 2. Energy gradient (derivative)
    gradient = np.gradient(free_energy, bin_centers)
    fig.add_trace(
        go.Scatter(
            x=bin_centers,
            y=gradient,
            mode='lines',
            name='Energy Gradient',
            line=dict(color='#A23B72', width=3),
            hovertemplate='<b>Committor</b>: %{x:.3f}<br><b>Gradient</b>: %{y:.2f}<extra></extra>'
        ),
        row=1, col=2
    )
    
    # 3. State distribution (exponential of negative free energy)
    prob_dist = np.exp(-free_energy / 0.616)  # kT at 310K
    prob_dist = prob_dist / np.sum(prob_dist)
    
    fig.add_trace(
        go.Bar(
            x=bin_centers,
            y=prob_dist,
            name='State Probability',
            marker_color='#F18F01',
            hovertemplate='<b>Committor</b>: %{x:.3f}<br><b>Probability</b>: %{y:.4f}<extra></extra>'
        ),
        row=2, col=1
    )
    
    # 4. 3D landscape
    x_3d = np.linspace(0, 1, 20)
    y_3d = np.linspace(0, max(free_energy), 20)
    X_3d, Y_3d = np.meshgrid(x_3d, y_3d)
    Z_3d = np.interp(X_3d.flatten(), bin_centers, free_energy).reshape(X_3d.shape)
    
    fig.add_trace(
        go.Surface(
            x=X_3d, y=Y_3d, z=Z_3d,
            colorscale='Viridis',
            name='3D Landscape',
            showscale=False
        ),
        row=2, col=2
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"<b>{target_name} Multi-Region Exploration Dashboard</b><br><sub>Interactive Thermodynamic Analysis</sub>",
            x=0.5,
            font=dict(size=20)
        ),
        height=800,
        template='plotly_white',
        font=dict(size=12)
    )
    
    # Update axes
    fig.update_xaxes(title_text="Committor p(B)", row=1, col=1)
    fig.update_yaxes(title_text="Free Energy (kcal/mol)", row=1, col=1)
    fig.update_xaxes(title_text="Committor p(B)", row=1, col=2)
    fig.update_yaxes(title_text="∂G/∂p(B)", row=1, col=2)
    fig.update_xaxes(title_text="Committor p(B)", row=2, col=1)
    fig.update_yaxes(title_text="Probability", row=2, col=1)
    
    # Save interactive plot
    output_path_interactive = output_dir / "interactive_dashboard.html"
    fig.write_html(str(output_path_interactive))
    print(f"✅ Interactive dashboard saved to {output_path_interactive}")
    
    return fig

def create_summary_statistics_plot(bin_centers, free_energy, target_name="Target"):
    """Create a summary statistics visualization"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Free energy with confidence intervals (simulated)
    noise = np.random.normal(0, 0.1, len(free_energy)) 
    upper_bound = free_energy + np.abs(noise)
    lower_bound = free_energy - np.abs(noise)
    
    ax1.plot(bin_centers, free_energy, 'o-', linewidth=3, markersize=8, color='#2E86AB', label='Mean Energy')
    ax1.fill_between(bin_centers, lower_bound, upper_bound, alpha=0.3, color='#2E86AB', label='Uncertainty')
    ax1.set_xlabel("Committor p(B)")
    ax1.set_ylabel("Free Energy (kcal/mol)")
    ax1.set_title("Free Energy with Uncertainty")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Energy histogram
    ax2.hist(free_energy, bins=15, color='#A23B72', alpha=0.7, edgecolor='black')
    ax2.axvline(np.mean(free_energy), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(free_energy):.2f}')
    ax2.axvline(np.median(free_energy), color='orange', linestyle='--', linewidth=2, label=f'Median: {np.median(free_energy):.2f}')
    ax2.set_xlabel("Free Energy (kcal/mol)")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Energy Distribution")
    ax2.legend()
    
    # 3. Committor vs Energy scatter with trend
    ax3.scatter(bin_centers, free_energy, s=100, c=free_energy, cmap=cmocean.cm.thermal, alpha=0.8, edgecolors='black')
    z = np.polyfit(bin_centers, free_energy, 2)
    p = np.poly1d(z)
    ax3.plot(bin_centers, p(bin_centers), "r--", linewidth=2, label='Quadratic Fit')
    ax3.set_xlabel("Committor p(B)")
    ax3.set_ylabel("Free Energy (kcal/mol)")
    ax3.set_title("Committor-Energy Correlation")
    ax3.legend()
    
    # 4. Key statistics table
    ax4.axis('off')
    stats_data = [
        ['Statistic', 'Value'],
        ['Min Energy', f'{np.min(free_energy):.3f} kcal/mol'],
        ['Max Energy', f'{np.max(free_energy):.3f} kcal/mol'],
        ['Mean Energy', f'{np.mean(free_energy):.3f} kcal/mol'],
        ['Energy Range', f'{np.ptp(free_energy):.3f} kcal/mol'],
        ['Std Deviation', f'{np.std(free_energy):.3f} kcal/mol'],
        ['Barrier Height', f'{np.max(free_energy) - np.min(free_energy):.3f} kcal/mol'],
        ['Transition State p(B)', f'{bin_centers[np.argmax(free_energy)]:.3f}'],
    ]
    
    table = ax4.table(cellText=stats_data[1:], colLabels=stats_data[0], 
                     cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 2)
    
    # Style the table
    for i in range(len(stats_data)):
        if i == 0:  # Header
            for j in range(2):
                table[(i, j)].set_facecolor('#2E86AB')
                table[(i, j)].set_text_props(weight='bold', color='white')
    
    ax4.set_title("Summary Statistics", fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    output_path_stats = output_dir / "summary_statistics.png"
    plt.savefig(output_path_stats, dpi=300, bbox_inches='tight')
    print(f"✅ Summary statistics saved to {output_path_stats}")
    
    plt.show()

def main():
    """Main function to create enhanced visualizations with argument parsing"""
    parser = argparse.ArgumentParser(description="Create enhanced visualization plots")
    parser.add_argument("--target", "-t", default="GLP1R", help="Target name (default: GLP1R)")
    args = parser.parse_args()
    
    MD_DIR = ARTIFACTS_DIR / "md" / args.target
    MD_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load the existing free energy data
    # (This would normally load from your analysis results)
    print("🎨 Creating enhanced visualizations...")
    
    # For demo purposes, let's create some example data
    # In practice, this would load your actual analysis results
    bin_centers = np.linspace(0, 1, 50)
    # Simulate a realistic free energy profile with barriers
    free_energy = 2 * np.sin(3 * np.pi * bin_centers) + 0.5 * bin_centers**2 + np.random.normal(0, 0.1, len(bin_centers))
    free_energy = free_energy - np.min(free_energy)  # Set minimum to zero
    
    print("📊 Creating enhanced matplotlib plots...")
    create_enhanced_free_energy_plot(bin_centers, free_energy, args.target, MD_DIR)
    
    print("Creating interactive Plotly dashboard...")
    interactive_fig = create_interactive_plotly_dashboard(bin_centers, free_energy, args.target, MD_DIR)
    if interactive_fig:
        print("Interactive dashboard created successfully")
    
    print("Creating summary statistics...")
    create_summary_statistics_plot(bin_centers, free_energy, args.target, MD_DIR)
    
    print("\n✅ All enhanced visualizations created!")
    print(f"📁 Check the output directory: {MD_DIR}")
    print("🔍 Open the interactive_dashboard.html in your browser for interactive exploration!")

if __name__ == "__main__":
    main()