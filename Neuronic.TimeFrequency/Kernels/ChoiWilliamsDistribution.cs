using System;
using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.TimeFrequency.Kernels
{
    public class ChoiWilliamsDistribution: DopplerLagKernel
    {
        public double Sigma { get; set; }

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