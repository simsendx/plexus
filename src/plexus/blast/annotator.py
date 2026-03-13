from collections import namedtuple

import pandas as pd


class BlastResultsAnnotator:
    def __init__(self, blast_df, target_map=None):
        """
        Annotate tabular results produced by running BLAST with
        `-outfmt 6`.

        params
            blast_df: DataFrame (n_hits, n_columns)
                Each row contains a BLAST hit.
            target_map: dict mapping synthetic primer IDs (e.g. "SEQ_3") to
                junction name(s) (e.g. "EXON1" or "EXON1|EXON2"). Optional.

        """
        self.blast_df = blast_df
        self.target_map = target_map or {}

    def build_annotation_dict(
        self,
        length_threshold=12,
        evalue_threshold=4,
        max_mismatches=2,
        three_prime_tolerance=0,
    ):
        """
        Build a dictionary specifying conditions for annotation.

        Parameters
        ----------
        three_prime_tolerance : int
            Maximum number of unaligned bases at the 3' end of the primer
            for a hit to still be considered "3'-anchored".  BLAST's local
            alignment clips terminal mismatches, so a primer with a mismatch
            at or near the 3' end will have ``qend < qlen``.  A tolerance
            of 2-3 catches these cases without accepting 5'-only hits.
        """
        self.annotation_params = {
            "length_threshold": length_threshold,
            "evalue_threshold": evalue_threshold,
            "max_mismatches": max_mismatches,
            "three_prime_tolerance": three_prime_tolerance,
        }
        # Keep annotation names for summarise_by_primer
        self.annotations = {
            "from_3prime": None,
            "length_pass_3prime": None,
            "evalue_pass_3prime": None,
            "predicted_bound": None,
        }

    def add_annotations(self):
        """
        Add annotations to `blast_df` using vectorized column operations.

        """
        p = self.annotation_params
        df = self.blast_df

        df["from_3prime"] = (df["qlen"] - df["qend"]) <= p["three_prime_tolerance"]
        df["length_pass_3prime"] = (
            (df["length"] >= p["length_threshold"])
            & (df["mismatch"] <= p["max_mismatches"])
            & df["from_3prime"]
        )
        df["evalue_pass_3prime"] = (df["evalue"] < p["evalue_threshold"]) & df[
            "from_3prime"
        ]
        df["predicted_bound"] = df["length_pass_3prime"] | df["evalue_pass_3prime"]

    def get_predicted_bound(self):
        """
        Get all BLAST matches where we predict that the query
        sequence will bind

        """
        if "predicted_bound" not in self.blast_df:
            raise ValueError(
                "`blast_df` must contain a column `predicted_bound`. Try running `.add_annotations()` first."
            )

        return self.blast_df.query("predicted_bound")

    def summarise_by_primer(self, output_path=None):
        """
        Summarise the BLAST results on a per-primer
        basis

        In essence, convert from a dataframe where each row
        is a BLAST hit, to a dataframe where each row
        is a primer, and columns give the total number of hits
        for that primer.

        """

        # Define structure to hold per-primer summaries
        PrimerBlastRecord = namedtuple(
            "PrimerBlastRecord",
            ["primer_name", "primer_pair_name", "target_name", "total_alignments"]
            + list(self.annotations),
        )

        # Collect summaries for every primer
        primer_blast_records = [
            PrimerBlastRecord(
                primer_name=qseqid,
                primer_pair_name=qseqid[:-2],
                target_name=self.target_map.get(qseqid, qseqid.split("_")[0]),
                total_alignments=qseqid_df.shape[0],
                **qseqid_df[list(self.annotations)].sum().to_dict(),
            )
            for qseqid, qseqid_df in self.blast_df.groupby("qseqid")
        ]

        # Convert to data frame
        self.blast_primer_df = pd.DataFrame(primer_blast_records).sort_values(
            "predicted_bound", ascending=False
        )

        # Optionally write
        if output_path is not None:
            self.blast_primer_df.to_csv(output_path, index=False)
