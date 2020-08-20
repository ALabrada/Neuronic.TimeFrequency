using System;

namespace Neuronic.TimeFrequency.Kernels
{
    public static class Windows
    {
        public static readonly WindowFunction Delta = new WindowFunction((i, n) => i == n/2 ? 1d : 0d); 
        
        public static readonly WindowFunction Rectangular = new WindowFunction((i, n) => 1d);

        public static readonly WindowFunction Bartlett = new WindowFunction((i, n) => i <= n/2 ? (2d*i)/n : (2d*(n - i))/n);

        public static readonly WindowFunction Hamming = new WindowFunction((i, n) => 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / n));

        public static readonly WindowFunction Hanning = new WindowFunction((i, n) => 0.5 * (1d - Math.Cos(2 * Math.PI * (i + 1) / (n + 2))));

        public static WindowFunction Gaussian(double alpha = 2.5) => new WindowFunction((i, n) =>
        {
            var e = 2 * alpha * (i - n / 2d) / n;
            return Math.Exp(-0.5 * e * e);
        });
    }
}