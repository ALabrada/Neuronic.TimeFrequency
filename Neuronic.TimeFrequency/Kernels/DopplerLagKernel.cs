using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.TimeFrequency.Kernels
{
    /// <summary>
    /// Base class for Doppler-lag kernels for time-frequency distributions.
    /// </summary>
    public abstract class DopplerLagKernel
    {
        /// <summary>
        /// Evaluates the kernel with the specified number of samples.
        /// </summary>
        /// <param name="count">The sample count.</param>
        /// <returns>The Doppler-lag smoothing kernel.</returns>
        public virtual double[,] Evaluate(int count)
        {
            var lagSample = 1;
            var dopplerSample = 1d / count;
            var g = new double[count, 2 * count];

            Evaluate(g, dopplerSample, lagSample);

            return g;
        }

        /// <summary>
        /// When implemented, evaluates the kernel in the specified buffer.
        /// </summary>
        /// <param name="g">The buffer.</param>
        /// <param name="dopplerSample">The doppler sample.</param>
        /// <param name="lagSample">The lag sample.</param>
        protected abstract void Evaluate(double[,] g, double dopplerSample, double lagSample);
    }
}