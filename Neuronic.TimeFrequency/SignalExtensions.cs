using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Neuronic.TimeFrequency
{
    static class SignalExtensions
    {
        public static Signal<Complex> ToComplex(this Signal<double> signal)
        {
            return new Signal<Complex>(signal.Samples.ToComplex(), signal.Delay, signal.SamplingRate);
        }

        public static void Integrate(this Signal<Complex> x)
        {
            var sum = Complex.Zero;
            for (int i = 0; i < x.Count; i++)
            {
                sum += x[i];
                x[i] = sum / (i * x.SamplingPeriod);
            }
        }
    }
}
