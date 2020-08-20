using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Choi-Williams Distribution.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Kernels.DopplerLagKernel" />
    public class ChoiWilliamsDistribution: DopplerLagKernel
    {
        /// <summary>
        /// Gets or sets the sigma parameter.
        /// </summary>
        public double Sigma { get; set; }

        /// <summary>
        /// Evaluates the kernel in the specified buffer.
        /// </summary>
        /// <param name="g">The buffer.</param>
        /// <param name="dopplerSample">The doppler sample.</param>
        /// <param name="lagSample">The lag sample.</param>
        protected override void Evaluate(double[,] g, double dopplerSample, double lagSample)
        {
            var halfRow = g.GetLength(0) / 2;
            var halfCol = (g.GetLength(1) - 1) / 2;

            for (int j = 0; j < g.GetLength(1); j++)
                g[0, j] = 1d;
            for (int i = 0; i < g.GetLength(0); i++)
                g[i, 0] = 1d;
            g[0, (g.GetLength(1) - 1) / 2 + 1] = 0d;

            var factor = 4 * Math.PI * Math.PI / Sigma;
            for (int m = 1; m <= halfCol; m++)
            for (int u = 1; u <= halfRow; u++)
            {
                var ul = u * dopplerSample;
                var ml = m * lagSample;
                var value = Math.Exp(-factor * (ul * ml) * (ul * ml));

                g[u, m] = value;
                g[g.GetLength(0) - u, m] = value;
                g[u, g.GetLength(1) - m] = value;
                g[g.GetLength(0) - u, g.GetLength(1) - m] = value;
            }
        }
    }
}