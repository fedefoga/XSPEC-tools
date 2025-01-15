from scipy.stats import beta
import numpy as np

def avevar(data):
    """
    Compute the average and variance of a dataset.
    Equivalent to `avevar` in the C++ code.
    """
    ave = np.mean(data)
    var = np.var(data, ddof=1)  # Use ddof=1 for sample variance
    return ave, var

def ftest(data1, data2):
    """
    Given the arrays data1[0..n1-1] and data2[0..n2-1], this routine returns 
    the value of f, and its p-value as prob. Small values of prob indicate that 
    the two arrays have significantly different variances.

    Parameters:
    data1, data2 : array-like
        Input datasets to compare.

    Returns:
    f : float
        The F-test statistic.
    prob : float
        The p-value associated with the F-test.
    """
    n1 = len(data1)
    n2 = len(data2)

    # Calculate average and variance for each dataset
    ave1, var1 = avevar(data1)
    ave2, var2 = avevar(data2)

    # Make F the ratio of the larger variance to the smaller one
    if var1 > var2:
        f = var1 / var2
        df1 = n1 - 1
        df2 = n2 - 1
    else:
        f = var2 / var1
        df1 = n2 - 1
        df2 = n1 - 1

    # Calculate probability using the incomplete beta function
    prob = 2.0 * beta.cdf(df2 / (df2 + df1 * f), 0.5 * df2, 0.5 * df1)
    if prob > 1.0:
        prob = 2.0 - prob

    return f, prob

