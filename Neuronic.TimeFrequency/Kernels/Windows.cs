using System;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Contains common window functions.
    /// </summary>
    public static class Windows
    {
        /// <summary>
        /// The Delta function is <c>1</c> at the center and <c>0</c> elsewhere.
        /// </summary>
        public static readonly WindowFunction Delta = new WindowFunction((i, n) => i == n/2 ? 1d : 0d);

        /// <summary>
        /// The Rectangular window.
        /// </summary>
        public static readonly WindowFunction Rectangular = new WindowFunction((i, n) => 1d);

        /// <summary>
        /// The Bartlett window.
        /// </summary>
        public static readonly WindowFunction Bartlett = new WindowFunction((i, n) => i <= n/2 ? (2d*i)/n : (2d*(n - i))/n);

        /// <summary>
        /// The Hamming window.
        /// </summary>
        public static readonly WindowFunction Hamming = new WindowFunction((i, n) => 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / n));

        /// <summary>
        /// The Hanning window.
        /// </summary>
        public static readonly WindowFunction Hanning = new WindowFunction((i, n) => 0.5 * (1d - Math.Cos(2 * Math.PI * (i + 1) / (n + 2))));

        /// <summary>
        /// Gets the Gaussian window with the specified alpha.
        /// </summary>
        /// <param name="alpha">The alpha.</param>
        /// <returns>The window.</returns>
        public static WindowFunction Gaussian(double alpha = 2.5) => new WindowFunction((i, n) =>
        {
            var e = 2 * alpha * (i - n / 2d) / n;
            return Math.Exp(-0.5 * e * e);
        });
    }
}