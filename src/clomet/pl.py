"""Diagnostic plots: the Leiden scan, choosing k, the splits and the consensus matrix."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram

no_date = {'CreationDate': None}



def resolution_scan(scan, used):
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.5))
    ax = axes[0, 0]
    ax.plot(scan.resolution, scan.silhouette, c='steelblue', lw=1, label='all resolutions')
    ax.scatter(used.resolution, used.silhouette, c='tomato', s=10, zorder=3, label='used (one per k)')
    ax.set_xlabel('Resolution'); ax.set_ylabel('Silhouette'); ax.legend()
    ax = axes[0, 1]
    ax.plot(scan.resolution, scan.k, c='steelblue', lw=1, marker='o', ms=2)
    ax.set_xlabel('Resolution'); ax.set_ylabel('Number of clusters (k)')
    ax = axes[1, 0]
    ax.plot(scan.resolution, scan.nmi_to_prev, c='darkorange', lw=1, marker='o', ms=2)
    ax.set_xlabel('Resolution'); ax.set_ylabel('NMI to previous resolution')
    ax = axes[1, 1]
    ax.scatter(used.k, used.silhouette, c='tomato', s=10)
    ax.set_xlabel('Number of clusters (k)'); ax.set_ylabel('Silhouette')
    fig.tight_layout()
    return fig


def choose_k(cuts, k_summary, k_range, chosen):
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.2))
    for k in k_range:
        sizes = cuts[f'k_{k}'].value_counts()
        axes[0].scatter([k] * len(sizes), sizes.values, s=6, c='steelblue', alpha=0.6, lw=0)
    axes[0].set_ylabel('Cluster size (clones)')
    axes[1].plot(k_range, k_summary.min_size, 'o-', ms=3, label='min')
    axes[1].plot(k_range, k_summary.max_size, 's-', ms=3, label='max')
    axes[1].plot(k_range, k_summary.mean_size, '^-', ms=3, label='mean')
    axes[1].fill_between(k_range, k_summary.min_size, k_summary.max_size, alpha=0.15)
    axes[1].set_ylabel('Cluster size'); axes[1].legend()
    axes[2].plot(k_range, k_summary.cv_size, 'o-', ms=3, c='purple')
    axes[2].set_ylabel('Coefficient of variation of sizes')
    for ax in axes:
        ax.set_xlabel('k (tree cut)'); ax.set_xticks(k_range)
        ax.axvline(chosen, c='grey', ls='--', lw=0.6)
    fig.tight_layout()
    return fig


def cut_ari(ari):
    fig, ax = plt.subplots(figsize=(3.6, 3))
    sns.heatmap(ari, annot=True, fmt='.2f', cmap='RdYlGn', vmin=-0.2, vmax=1, ax=ax, annot_kws=dict(fontsize=5),
                cbar_kws=dict(label='Adjusted Rand index'))
    ax.set_xlabel('k'); ax.set_ylabel('k'); ax.set_title('Agreement between tree cuts')
    return fig


def cut_clustermaps(coassoc, tree, cuts, ks, pdf_path):
    """Consensus matrix with the clusters of each cut in `ks` as colour bars, one page per k."""
    with PdfPages(pdf_path, metadata=no_date) as pdf:
        for k in ks:
            lab = cuts[f'k_{k}']
            pal = dict(zip(sorted(lab.unique()), sns.color_palette('tab20', lab.nunique())))
            g = sns.clustermap(coassoc, row_linkage=tree, col_linkage=tree, row_colors=lab.map(pal),
                               col_colors=lab.map(pal), cmap='Blues', xticklabels=False, yticklabels=False,
                               figsize=(4.5, 4.5), cbar_kws=dict(label='P(same Leiden cluster)'))
            g.ax_heatmap.collections[0].set_rasterized(True)
            g.ax_row_colors.set_axis_off(); g.ax_col_colors.set_axis_off()
            g.ax_heatmap.set_xlabel(''); g.ax_heatmap.set_ylabel('')
            g.fig.suptitle(f'k = {k}', y=1.02)
            pdf.savefig(g.fig)
            plt.close(g.fig)


def sankey(cluster_assignments, k_range):
    import plotly.graph_objects as go
    nodes = []
    node_labels = []
    node_colors = []
    node_positions_x = []
    node_positions_y = []

    links_source = []
    links_target = []
    links_value = []
    links_color = []

    max_clusters = cluster_assignments.max().max()
    colors = sns.color_palette("tab20", max_clusters).as_hex()

    node_idx = 0
    node_map = {}  # (k, cluster_id) -> node_idx

    for i, k in enumerate(k_range):
        col_name = f"k_{k}"
        unique_clusters = sorted(cluster_assignments[col_name].unique())

        x_pos = i / (len(k_range) - 1) if len(k_range) > 1 else 0.5

        for j, cluster_id in enumerate(unique_clusters):
            size = (cluster_assignments[col_name] == cluster_id).sum()

            node_map[(k, cluster_id)] = node_idx
            nodes.append(node_idx)
            node_labels.append(f"K={k}<br>C{cluster_id}<br>n={size}")
            node_colors.append(colors[(cluster_id - 1) % len(colors)])
            node_positions_x.append(x_pos)
            node_positions_y.append(j / max(len(unique_clusters) - 1, 1))

            node_idx += 1

    k_list = list(k_range)
    for i in range(len(k_list) - 1):
        k_current = k_list[i]
        k_next = k_list[i + 1]

        col_current = f"k_{k_current}"
        col_next = f"k_{k_next}"

        flow_df = cluster_assignments[[col_current, col_next]].copy()
        flow_matrix = pd.crosstab(flow_df[col_current], flow_df[col_next])

        for cluster_current in flow_matrix.index:
            for cluster_next in flow_matrix.columns:
                value = flow_matrix.loc[cluster_current, cluster_next]
                if value > 0:
                    source_idx = node_map[(k_current, cluster_current)]
                    target_idx = node_map[(k_next, cluster_next)]

                    links_source.append(source_idx)
                    links_target.append(target_idx)
                    links_value.append(value)

                    hex_color = node_colors[source_idx]
                    r = int(hex_color[1:3], 16)
                    g = int(hex_color[3:5], 16)
                    b = int(hex_color[5:7], 16)
                    rgba_color = f'rgba({r},{g},{b},0.4)'
                    links_color.append(rgba_color)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels,
            color=node_colors,
            x=node_positions_x,
            y=node_positions_y
        ),
        link=dict(
            source=links_source,
            target=links_target,
            value=links_value,
            color=links_color
        )
    )])

    fig.update_layout(
        title="Cluster Flow Across Different K Values",
        font=dict(size=10),
        height=800,
        width=1200
    )

    return fig


def split(cid, sp, renumber, n_groups=None):
    best, tried, z, members, sim = sp.best, sp.tried, sp.z, sp.members, sp.sim
    side = {**{c: 'A' for c in best['a']}, **{c: 'B' for c in best['b']}}
    name = {s: renumber[f'{cid}{s}'] for s in 'AB'}
    col = {'A': '#2196F3', 'B': '#FF5722'}
    fig, axes = plt.subplots(1, 3, figsize=(10, 3))
    dn = dendrogram(z, ax=axes[0], no_labels=True, color_threshold=0, above_threshold_color='black',
                    link_color_func=lambda _: 'black')
    leaves = [members[i] for i in dn['leaves']]
    top = z[:, 2].max()
    axes[0].bar(np.arange(len(leaves)) * 10 + 5, top * 0.04, bottom=-top * 0.05, width=10, lw=0,
                color=[col[side[c]] for c in leaves])
    axes[0].set_ylim(-top * 0.05, top * 1.02)
    axes[0].legend(handles=[Patch(fc=col[s], label=name[s]) for s in 'AB'], title='Meta-clone', loc='upper right')
    group = f'Group {cid}' + (f' of {n_groups}' if n_groups else '')
    axes[0].set_title(f'{group}: tree of consensus distances (k = {best["k"]})')
    axes[0].set_ylabel('Distance')
    order = [c for c in leaves if side[c] == 'A'] + [c for c in leaves if side[c] == 'B']
    im = axes[1].imshow(sim.loc[order, order].values, cmap='Blues', aspect='auto', interpolation='nearest',
                        rasterized=True)
    n_a = len(best['a'])
    axes[1].axvline(n_a - 0.5, c='red', lw=0.8); axes[1].axhline(n_a - 0.5, c='red', lw=0.8)
    axes[1].set_xticks([]); axes[1].set_yticks([])
    axes[1].set_title(f'{name["A"]} (n={n_a}) | {name["B"]} (n={len(best["b"])})')
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04).set_label('Consensus similarity')
    s = axes[2].scatter([t['balance'] for t in tried], [t['separation'] for t in tried],
                        c=[t['score'] for t in tried], cmap='viridis', s=6)
    axes[2].scatter(best['balance'], best['separation'], s=60, c='red', marker='*', label='best')
    axes[2].axvline(0.15, c='green', ls='--', lw=0.8, label='min balance')
    axes[2].set_xlabel('Balance (B / A)'); axes[2].set_ylabel('Separation'); axes[2].legend()
    fig.colorbar(s, ax=axes[2]).set_label('Score')
    fig.tight_layout()
    return fig


def coassociation(coassoc, tree, meta_clone, colors):
    """Consensus matrix ordered by the tree, with the meta-clones as colour bars."""
    c = meta_clone.map(colors)
    g = sns.clustermap(coassoc, row_linkage=tree, col_linkage=tree, row_colors=c, col_colors=c, cmap='Blues',
                       xticklabels=False, yticklabels=False, figsize=(6, 6),
                       cbar_kws=dict(label='P(same Leiden cluster)'))
    g.ax_heatmap.collections[0].set_rasterized(True)
    g.ax_row_colors.set_axis_off(); g.ax_col_colors.set_axis_off()
    g.ax_heatmap.set_xlabel(''); g.ax_heatmap.set_ylabel('')
    g.ax_row_dendrogram.legend(handles=[Patch(fc=colors[m], label=m) for m in colors],
                               title='Meta-clone', bbox_to_anchor=(5.6, 1.1), loc='upper left', frameon=False)
    return g.fig


create_cluster_flow_sankey = sankey
