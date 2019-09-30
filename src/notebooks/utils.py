""" Utility functions, some of which are borrowed from Pat Walters: https://github.com/PatWalters/metk """

def rmse(pred_array, ref_array):
    """
    Calculate root mean squared (rms) error
    :param pred_array: the predicted values
    :param ref_array: the reference values
    :return: the rms error
    """
    return np.sqrt(np.mean((pred_array - ref_array) ** 2))


def mean_absolute_error(pred_array, ref_array):
    """
    Calculate mean absolute error
    :param pred_array: the predicted values
    :param ref_array: the reference values
    :return: the mean absolute error
    """
    return np.mean(np.abs(pred_array - ref_array))

def max_possible_correlation(vals, error=1 / 3.0, method=pearsonr, cycles=1000):
    """
    Calculate the maximum possible correlation given a particular experimental error
    Based on Brown, Muchmore, Hajduk http://www.sciencedirect.com/science/article/pii/S1359644609000403
    :param vals: experimental values (should be on a log scale)
    :param error: experimental error
    :param method: method for calculating the correlation, must take 2 lists and return correlation and p_value
    :param cycles: number of random cycles
    :return: maximum possible correlation
    """
    cor_list = []
    for i in range(0, cycles):
        noisy_vals = []
        for val in vals:
            noisy_vals.append(val + np.random.normal(0, error))
        cor_list.append(method(vals, noisy_vals)[0])
    return np.mean(cor_list)
