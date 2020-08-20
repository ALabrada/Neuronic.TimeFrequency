using System;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Smoothed Wigner-Ville Lag-Independent kernel.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Kernels.DopplerKernel" />
    public class SmoothedWignerVilleDistribution : DopplerKernel
    {
        /// <summary>
        /// Evaluates the kernel in the specified buffer.
        /// </summary>
        /// <param name="g">The buffer.</param>
        /// <param name="dopplerSample">The doppler sample.</param>
        /// <param name="lagSample">The lag sample.</param>
        protected override void Evaluate(double[,] g, double dopplerSample, double lagSample)
        {
            var g1 = new double[g.GetLength(0)];
            var count = Evaluate(g.GetLength(0), g.GetLength(1), g1, 0);
            if (UseDopplerDomain)
                PadWindow(new ArraySegment<double>(g1, 0, count), g1);

            for (int i = 0; i < g.GetLength(0); i++)
            for (int j = 0; j < g.GetLength(1); j++)
                g[i, j] = g1[i];
        }
    }
}