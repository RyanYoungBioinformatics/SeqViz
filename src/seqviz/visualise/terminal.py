from rich.console import Console
from rich.text import Text

console = Console()

def print_alignment(seq1_aligned, seq2_aligned, label: str | None = None, score: int | None = None):
    if label:
        console.rule(f"[bold cyan]{label}")

    if score is not None:
        console.print(f"[bold]Score:[/] {score}")

    line1 = Text()
    line2 = Text()
    for a, b in zip(seq1_aligned, seq2_aligned):
        if a == '-' or b == '-':
            line1.append(a, style="bold yellow")
            line2.append(b, style="bold yellow")
        elif a == b:
            line1.append(a, style="bold green")
            line2.append(b, style="bold green")
        else:
            line1.append(a, style="bold red")
            line2.append(b, style="bold red")
    console.print(line1)
    console.print(line2)
