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

        public static IEnumerable<T> Sample<T>(this Signal<T> psi, double freq)
        {
            for (int i = 0; i < freq * psi.Count * psi.SamplingPeriod; i++)
            {
                var j = (int)Math.Floor(i / (freq * psi.SamplingPeriod));
                yield return psi[j];
            }
        }

        public static void Convolve(this IList<float> samples, IList<Complex> kernel, IList<Complex> result)
        {
            for (int i = 0; i < result.Count; i++)
            {
                var sum = Complex.Zero;
                var start = i - (result.Count + 1) / 2 - (kernel.Count + 1) / 2 + (samples.Count + 1) / 2;
                for (var j = Math.Max(0, -start); j < kernel.Count && j + start < samples.Count; j++)                
                    sum += samples[j + start] * kernel[kernel.Count - 1 - j];
                result[i] = sum;
            }
        }
    }
}
