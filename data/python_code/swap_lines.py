#!/usr/bin/env python
import argparse
import random

def main():
    parser = argparse.ArgumentParser(
        description="Randomly swap (shuffle) the lines of a file."
    )
    parser.add_argument("input_file", type=str, help="Path to the input file.")
    parser.add_argument("output_file", type=str, help="Path to the output file.")
    parser.add_argument(
        "--swap-only",
        action="store_true",
        help="If set, only swap two random lines instead of shuffling all lines.",
    )
    args = parser.parse_args()

    # Read the file's lines.
    with open(args.input_file, "r") as infile:
        lines = infile.readlines()

    # Perform swapping.
    if args.swap_only:
        if len(lines) < 2:
            print("Not enough lines to swap.")
        else:
            i, j = random.sample(range(len(lines)), 2)
            lines[i], lines[j] = lines[j], lines[i]
    else:
        random.shuffle(lines)

    # Write the swapped lines to the output file.
    with open(args.output_file, "w") as outfile:
        outfile.writelines(lines)
    print(f"Processed {len(lines)} lines. Output written to {args.output_file}")

if __name__ == "__main__":
    main()
