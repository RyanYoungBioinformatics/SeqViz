import seaborn as sns
import matplotlib.pyplot as plt

def plot_matrix(matrix, seq1, seq2, output_path=None):
    ax = sns.heatmap(
        matrix,
        xticklabels=['-'] + list(seq2),
        yticklabels=['-'] + list(seq1),
        annot=True,      # shows score values in each cell
        fmt='d',
        cmap='coolwarm'
    )
    plt.title("Alignment Scoring Matrix")
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    else:
        plt.show()
