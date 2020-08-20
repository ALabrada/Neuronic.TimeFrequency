namespace Neuronic.TimeFrequency.Kernels
{
    public class WignerVilleDistribution : DopplerLagKernel
    {
        protected override void Evaluate(double[,] g, double dopplerSample, double lagSample)
        {
            for (int i = 0; i < g.GetLength(0); i++)
            for (int j = 0; j < g.GetLength(1); j++)
                g[i, j] = 1d;
        }
    }
}