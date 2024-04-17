#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    def __init__(
        self,
        individual_col="individual",
        sample_col="sample",
        modbam_col="modbam_5mCG",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            individual_col (str): The name of the column that contains the individual name
                (default "individual").
            sample_col (str): The name of the column that specifies the colon-delimited list of tissues of the sample.
            modbam_col (str): The modified BAM file associated with the patient and sample.
        """
        super().__init__(**kwargs)
        self._individual_col = individual_col
        self._sample_col = sample_col
        self._modbam_col = modbam_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_individual(row)
        self._validate_sample(row)
        self._validate_modbam(row)
        self._seen.add((row[self._individual_col], row[self._sample_col], row[self._modbam_col]))
        self.modified.append(row)

    def _validate_individual(self, row):
        """Assert that the individual name exists and convert spaces to underscores."""
        if len(row[self._individual_col]) <= 0:
            raise AssertionError("Individual input is required.")
        # Sanitize individuals slightly.
        row[self._individual_col] = row[self._individual_col].replace(" ", "_")

    def _validate_sample(self, row):
        """Assert that the tissue names exists."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("The samples column is required.")
    
    def _validate_modbam(self, row):
        """Assert that the column with modified bam files exists."""
        if len(row[self._modbam_col]) <= 0:
            raise AssertionError("Modified BAMs is required.")

    def validate_unique_individuals(self):
        """
        Assert that the combination of individual name and modified BAM file is unique.

        In addition to the validation, also rename all individuals to have a suffix of _T{n}, where n is the
        number of times the same individual exist, but with different folders, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of individual name and modified BAM directory must be unique.")
        seen = Counter()
        for row in self.modified:
            individual = row[self._individual_col]
            seen[individual] += 1
            row[self._individual_col] = f"{individual}_T{seen[individual]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,raw_dir
            SAMPLE_1,run_dir_1
            SAMPLE_2,run_dir_2

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    required_columns = {"individual", "sample", "modbam_5mCG"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        # checker.validate_unique_individuals() # NB this may be needed later
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
