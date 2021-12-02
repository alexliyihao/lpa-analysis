import pandas as pd
import seaborn as sns

def encode(df):
    """
    The function one-hot encoded the DataFrame
    Args:
        df: pandas.DataFrame, the output from read_vcf
    Return:
        df: pandas.DataFrame, the DataFrame which is one-hot-encoded
    """
    df = pd.get_dummies(df, drop_first = True, prefix_sep = "-")
    # from alphanumerical order all the 0/0 should be dropped, just in case
    assert(len([i for i in df.columns if "0/0" in i]) == 0)
    return df


def edav(snp_ys,label):
    print(f"----------------Summary for label {label}----------------")
    for snp_y in snp_ys:
        print(f"NA exist:{snp_y[label].isna().any()}")
        vc = snp_y[label].value_counts()
        if len(vc) > 20:
            sns.displot(snp_y[label])
        else:
            print(vc)
