import scipy.integrate as integrate
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt
import time
import multiprocessing
import numpy as np


class CustomRvContinuous(rv_continuous):
    """
    Custom Continuous distribution based on probability distribution function string
    """
    def __init__(self, pdf_str, min_value, max_value):
        """
        Constructs custom continuous random variable
        @param pdf_str: Function string with "x" as a function parameter. Can use np for numpy
        @param min_value: Lower bound
        @param max_value: Upper bound
        """
        assert min_value < max_value, "min value must be < max value but instead min_value = {}, max_value = {}".format(
            min_value, max_value)
        super().__init__()
        self._min_value = min_value
        self._max_value = max_value
        self._raw_pdf = eval("lambda x: {func_str}".format(func_str=pdf_str))
        self._integral_value, self._integral_error = integrate.quad(self._raw_pdf, self._min_value, self._max_value)

    # We should implement cumulative distribution function, because the generic implementation is very slow
    def _cdf(self, x, *args):
        if x <= self._min_value:
            return 0
        if x >= self._max_value:
            return 1

        value, error = integrate.quad(self._raw_pdf, self._min_value, x)
        return value / self._integral_value

    def _pdf(self, x, **kwargs):
        # If value outside the boundaries pdf must return 0
        if x < self._min_value or x > self._max_value:
            return 0

        # Otherwise, normalize raw_pdf with the integral value
        return self._raw_pdf(x) / self._integral_value

    @staticmethod
    def create_portion_sample_rvs(pdf_str, min_value, max_value, portion_size):
        """
        Single job of probability distribution sampling.
        @note
        Each job should create a separate custom rv class and random state to ensure we don't have any locks that will
        reduce efficiency.
        @param pdf_str: Probability distribution function string
        @param min_value: Lower bound
        @param max_value: Upper bound
        @param portion_size: Sample size
        @return: np.ndarray array of size <portion_size> with values sampled based on the probability distribution
        """
        custom_rv = CustomRvContinuous(pdf_str=pdf_str, min_value=min_value, max_value=max_value)
        return custom_rv.rvs(size=portion_size, random_state=np.random.RandomState())

    @staticmethod
    def multithreaded_rvs(size, pdf_str, min_value, max_value):
        """
        Multi threaded sampling from a certain probability distribution
        @param size: Sample size
        @param pdf_str: Probability distribution function string
        @param min_value: Lower bound
        @param max_value: Upper bound
        @return: np.ndarray array of size <size> with values sampled based on the probability distribution
        """
        cpu_count = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=cpu_count)
        async_results = []

        for i in range(cpu_count - 1):
            async_results.append(pool.apply_async(CustomRvContinuous.create_portion_sample_rvs, (pdf_str, min_value, max_value, int(size / cpu_count))))
        async_results.append(pool.apply_async(CustomRvContinuous.create_portion_sample_rvs, (pdf_str, min_value, max_value, size - (cpu_count - 1) * int(size / cpu_count))))

        samples = async_results[0].get()
        for i in range(1, cpu_count):
            cur_samples = async_results[i].get()
            samples = np.concatenate((samples, cur_samples))
        return samples


if __name__ == '__main__':
    FUNC_STRING = "x ** 2"
    NORM_FUNC_STRING = "np.exp(-x**2 / 2.) / np.sqrt(2.0 * np.pi)"

    start_time = time.time()
    arr = CustomRvContinuous.multithreaded_rvs(size=500000, pdf_str=NORM_FUNC_STRING, min_value=-5, max_value=5)
    print("Took {}".format(time.time() - start_time))
    plt.hist(arr, 50)
    plt.show()
