using System;

namespace Neuronic.TimeFrequency.Kernels
{
    public class SmoothedWignerVilleDistribution : DopplerKernel
    {
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