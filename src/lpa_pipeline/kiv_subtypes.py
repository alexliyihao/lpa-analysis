"""
The KIV-2 types A, B and C are defined by three synonymous SNPs
at positions 14, 41 and 86 of the first KIV-2 exon which can be seen via
https://www.atherosclerosis-journal.com/cms/10.1016/j.atherosclerosis.2022.04.003/attachment/32899cc9-f3a8-4354-86a8-5ebcec782688/mmc1.pdf

The file from pipeline: https://github.com/genepi/lpa-pipeline/blob/master/files/maplocus_v3.txt
suggested that Exon 421 should be bp 581-740 from his pipeline. Which can be confirmed by his next sentence
"A T/C SNP 119 bp downstream of exon 1 (position 859 in [1])" Thus pos 14(A or G), 41(T or C), 86(A or T)
should be pos 594, 621 and 666 in our pipeline, respectively.

Example::

    krs = KIVRepeatSubtype()
    krs.predict_subtype(encoding_result = encoding_result)

where encoding result is a pd.DataFrame with all value binary(or equivalent)
each row should stand for an individual(id as index) and each column stand for a variant,
using "pos-ref/alt" format. all of ["594-A/G", "621-T/C", "666-A/T"] have to be in the columns
"""
import pandas as pd


class KIVRepeatSubtype:

    def __init__(self) -> None:
        self.encoding_result = None
        self.related_SNPs = ["594-A/G", "621-T/C", "666-A/T"]

    def predict_subtype(self, encoding_result: pd.DataFrame) -> pd.DataFrame:
        """predict the KIV‐2 types A, B and C, wrappers for _predict_subtype

        Args:
            encoding_result: pd.DataFrame with all value binary or equivalent
            each row should stand for an individual, where id is used as index
            and each column stand for a variant, using "pos-ref/alt" format
        Returns:
            pd.DataFrame: have the same index as encoding_result and subtype prediction
            using column name "KIV-2_subtypes" as the only column
        """
        assert encoding_result.columns.isin(self.related_SNPs).sum() == 3, \
            "the encoding result doesn't have all 3 SNPs necessary"
        self.encoding_result = encoding_result
        encoding_result_subtype = self.encoding_result.loc[
                                  :, encoding_result.columns.isin(self.related_SNPs)
                                  ].apply(self._predict_subtype, axis=1)
        print(encoding_result_subtype.value_counts())
        return pd.DataFrame(encoding_result_subtype, columns=["KIV-2_subtypes"])

    def _predict_subtype(self, x: pd.DataFrame) -> str:
        """predict the KIV‐2 types A, B and C described by:

        https://doi.org/10.1016/j.atherosclerosis.2022.04.003
        """
        if (x["594-A/G"] == 0) & (x["621-T/C"] == 0) & (x["666-A/T"] == 0):
            return "A"
        elif (x["594-A/G"] == 1) & (x["621-T/C"] == 1) & (x["666-A/T"] == 1):
            return "B"
        elif (x["594-A/G"] == 1) & (x["621-T/C"] == 1) & (x["666-A/T"] == 0):
            return "C"
        else:
            return "unspecified"
