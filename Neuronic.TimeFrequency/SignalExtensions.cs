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
                x[i] = sum * x.SamplingPeriod;
            }
        }

        public static void Conjugate(this Signal<Complex> x)
        {
            for (int i = 0; i < x.Count; i++)
                x[i] = Complex.Conjugate(x[i]);
        }

        public static void Differentiate(this Signal<Complex> x)
        {
            for (int i = 1; i < x.Count; i++)
                x[i - 1] = x[i] - x[i - 1];
        }

        public static void Convolve(this float[] samples, Complex[] kernel, Complex[] result)
        {
            for (int i = 0; i < result.Length; i++)
            {
                var sum = Complex.Zero;
                var start = i - (result.Length + 1) / 2 - (kernel.Length + 1) / 2 + (samples.Length + 1) / 2;
                for (var j = Math.Max(0, -start); j < kernel.Length && j + start < samples.Length; j++)                
                    sum += samples[j + start] * kernel[kernel.Length - 1 - j];
                result[i] = sum;
            }
        }
    }
}
