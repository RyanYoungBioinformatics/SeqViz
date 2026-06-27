import click
from seqviz.scoring import MATCH, MISMATCH, GAP
from seqviz.algorithms import needleman_wunsch as nw
from seqviz.algorithms import smith_waterman as sw
from seqviz.fasta_utils import read_fasta
from seqviz.visualise.heatmap import plot_matrix
from seqviz.visualise.terminal import print_alignment


@click.command()
@click.option(
    "--seq1",
    default=None,
    help="First DNA or protein sequence e.g. GATTACA (or use --file1)",
)
@click.option(
    "--seq2",
    default=None,
    help="Second DNA or protein sequence e.g. GCATGCU (or use --file2)",
)
@click.option(
    "--file1",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to a FASTA file containing the first sequence",
)
@click.option(
    "--file2",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to a FASTA file containing the second sequence",
)
@click.option(
    "--algorithm",
    default="nw",
    type=click.Choice(["nw", "sw", "both"]),
    show_default=True,
    help="nw = global alignment, sw = local alignment, both = run both",
)
@click.option(
    "--match",
    "match_score",        # python variable name — avoids clash with built-in
    default=MATCH,
    type=int,
    show_default=True,
    help="Score awarded for a matching character",
)
@click.option(
    "--mismatch",
    "mismatch_score",
    default=MISMATCH,
    type=int,
    show_default=True,
    help="Penalty for a mismatching character",
)
@click.option(
    "--gap",
    "gap_penalty",
    default=GAP,
    type=int,
    show_default=True,
    help="Gap penalty applied when a gap is introduced",
)
@click.option(
    "--output",
    default=None,
    type=click.Path(),
    help="Optional file path to save the scoring matrix heatmap  e.g. matrix.png",
)
@click.option(
    "--mode",
    default="dna",
    type=click.Choice(["dna", "protein"]),
    show_default=True,
    help="Sequence mode — protein mode uses simple scoring for now",
)
def main(seq1, seq2, file1, file2, algorithm, match_score, mismatch_score, gap_penalty, output, mode):
    """
    SeqViz — Sequence Alignment Visualiser

    Align two DNA or protein sequences using Needleman-Wunsch (global)
    and/or Smith-Waterman (local) algorithms.

    Inputs may be provided as inline sequences via --seq1/--seq2 or as FASTA files
    via --file1/--file2.
    """
    if file1:
        seq1 = read_fasta(file1)
    elif not seq1:
        raise click.UsageError("Provide --seq1 or --file1")
    else:
        seq1 = seq1.upper()

    if file2:
        seq2 = read_fasta(file2)
    elif not seq2:
        raise click.UsageError("Provide --seq2 or --file2")
    else:
        seq2 = seq2.upper()

    if mode == "protein":
        click.echo(
            "\nNote: protein mode currently uses simple match/mismatch scoring. "
            "BLOSUM62 support coming soon."
        )

    # Collect results so we can pass them to visualisation later
    results = {}

    if algorithm in ("nw", "both"):
        a1, a2, matrix, score = nw.align(
            seq1, seq2,
            match=match_score,
            mismatch=mismatch_score,
            gap=gap_penalty,
        )
        results["Needleman-Wunsch (Global)"] = (a1, a2, matrix, score)

    if algorithm in ("sw", "both"):
        a1, a2, matrix, score = sw.align(
            seq1, seq2,
            match=match_score,
            mismatch=mismatch_score,
            gap=gap_penalty,
        )
        results["Smith-Waterman (Local)"] = (a1, a2, matrix, score)

    # Print each result
    #for label, (a1, a2, matrix, score) in results.items():
        #print_alignment(label, a1, a2, score)

    # Heatmap stub
    for label, (a1, a2, matrix, score) in results.items():
        print_alignment(a1, a2, label=label, score=score)

        if output:
            if algorithm == "both":
                suffix = "_nw" if "Needleman" in label else "_sw"
                base, ext = output.rsplit(".", 1)
                save_path = f"{base}{suffix}.{ext}"
            else:
                save_path = output

            saved_to = plot_matrix(matrix, seq1, seq2, title=label, output_path=save_path)
            click.echo(f"\n  Heatmap saved → {saved_to}")
        else:
            plot_matrix(matrix, seq1, seq2, title=label)


if __name__ == "__main__":
    main()
