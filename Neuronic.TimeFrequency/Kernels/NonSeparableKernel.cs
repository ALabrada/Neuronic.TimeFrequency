using System.Collections.Generic;
using System.Numerics;

namespace Neuronic.TimeFrequency.Kernels
{
    public abstract class NonSeparableKernel
    {
        public virtual double[,] Evaluate(int count)
        {
            var lagSample = 1;
            var dopplerSample = 1d / count;
            var g = new double[count, 2 * count];

            Evaluate(g, dopplerSample, lagSample);

            return g;
        }

        protected abstract void Evaluate(double[,] g, double dopplerSample, double lagSample);
    }
}