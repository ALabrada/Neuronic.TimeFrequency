using System;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Product kernel
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Kernels.DopplerKernel" />
    public class ProductKernel : DopplerKernel
    {
        /// <summary>
        /// Evaluates the kernel in the specified buffer.
        /// </summary>
        /// <param name="g">The buffer.</param>
        /// <param name="dopplerSample">The doppler sample.</param>
        /// <param name="lagSample">The lag sample.</param>
        protected override void Evaluate(double[,] g, double dopplerSample, double lagSample)
        {
            var lDop = g.GetLength(0) * g.GetLength(1);
            var g1 = new double[lDop];
            var count = Evaluate(lDop, g.GetLength(1), g1, 0);
            if (UseDopplerDomain)
                PadWindow(new ArraySegment<double>(g1, 0, count), new ArraySegment<double>(g1, 0, g.GetLength(1)));

            for (int i = 0; i < g.GetLength(1); i++)
                g[0, i] = 1d;
            for (int i = 0; i < g.GetLength(0); i++)
                g[i, 0] = 1d;

            var lHalf = (g.GetLength(1) - 1) / 2;
            g[0, lHalf + 1] = 0d;

            for (int u = 1; u <= g.GetLength(0)/2; u++)
            for (int m = 1; m <= lHalf; m++)
            {
                var im = (u * m) % lDop;

                var value = g1[im];
                g[u, m] = value;
                g[g.GetLength(0) - u, m] = value;
                g[u, g.GetLength(1) - m] = value;
                g[g.GetLength(0) - u, g.GetLength(1) - m] = value;
            }
        }
    }
}