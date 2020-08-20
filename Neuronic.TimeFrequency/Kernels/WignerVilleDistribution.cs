namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// The Wigner-Ville Distribution.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Kernels.DopplerLagKernel" />
    public class WignerVilleDistribution : DopplerLagKernel
    {
        /// <summary>
        /// Evaluates the kernel in the specified buffer.
        /// </summary>
        /// <param name="g">The buffer.</param>
        /// <param name="dopplerSample">The doppler sample.</param>
        /// <param name="lagSample">The lag sample.</param>
        protected override void Evaluate(double[,] g, double dopplerSample, double lagSample)
        {
            for (int i = 0; i < g.GetLength(0); i++)
            for (int j = 0; j < g.GetLength(1); j++)
                g[i, j] = 1d;
        }
    }
}