using Accord.Math;
using System;
using System.Collections;
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
            return new Signal<Complex>(signal.Samples.ToComplex(), signal.Start, signal.SamplingRate);
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

        public static void Integrate(this Signal<double> x)
        {
            var sum = 0d;
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

        public static void Differentiate(this Signal<double> x)
        {
            for (int i = 1; i < x.Count; i++)
                x[i - 1] = x[i] - x[i - 1];
        }

        public static IEnumerable<T> Sample<T>(this IReadOnlySignal<T> psi, double freq)
        {
            for (int i = 0; i < freq * psi.Count * psi.SamplingPeriod; i++)
            {
                var j = (int)Math.Floor(i / (freq * psi.SamplingPeriod));
                yield return psi[j];
            }
        }

        public static void Convolve(this IReadOnlyList<double> samples, IReadOnlyList<Complex> kernel, IList<Complex> result, int stride = 1)
        {
            var dif = -(stride * result.Count) / 2 - (kernel.Count) / 2 + (samples.Count) / 2 + Math.Max(0, (samples.Count & 1) - (kernel.Count & 1));
            for (int i = 0; i < result.Count; i++)
            {
                var sum = Complex.Zero;
                var start = stride * i + dif;
                for (var j = Math.Max(0, -start); j < kernel.Count && j + start < samples.Count; j++)                
                    sum += samples[j + start] * kernel[kernel.Count - 1 - j];
                result[i] = sum;
            }
        }

        public static void Convolve(this IReadOnlyList<double> samples, IReadOnlyList<double> kernel, IList<double> result, int stride = 1)
        {
            var dif = -(stride * result.Count) / 2 - (kernel.Count) / 2 + (samples.Count) / 2 + Math.Max(0, (samples.Count & 1) - (kernel.Count & 1));
            for (int i = 0; i < result.Count; i++)
            {
                var sum = 0d;
                var start = stride * i + dif;
                for (var j = Math.Max(0, -start); j < kernel.Count && j + start < samples.Count; j++)
                    sum += samples[j + start] * kernel[kernel.Count - 1 - j];
                result[i] = sum;
            }
        }

        public static void CopyTo<T>(this IEnumerable<T> samples, T[] array, int start)
        {
            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (array == null) throw new ArgumentNullException(nameof(array));
            if (start < 0 || start >= array.Length) throw new ArgumentOutOfRangeException(nameof(start));
            foreach (var item in samples.Take(array.Length - start))
                array[start++] = item;
        }

        public static int CountOscilations(this IReadOnlyList<double> x)
        {
            var count = 0;
            var isTop = false;
            var isDown = false;
            for (int i = 1; i < x.Count - 1; i++)
            {
                if (x[i - 1] < x[i] && x[i + 1] < x[i])
                    count++;
                if (x[i - 1] > x[i] && x[i + 1] > x[i])
                    count++;

                if (x[i - 1] < x[i] && x[i + 1] == x[i])
                {
                    isTop = true;
                    isDown = false;
                }
                if (x[i - 1] == x[i] && x[i + 1] < x[i])
                {
                    if (isTop)
                        count++;
                    isTop = false;
                }

                if (x[i - 1] > x[i] && x[i + 1] == x[i])
                {
                    isTop = false;
                    isDown = true;
                }
                if (x[i - 1] == x[i] && x[i + 1] < x[i])
                {
                    if (isDown)
                        count++;
                    isDown = false;
                }
            }

            return count;
        }

        #region Casting
        public static IReadOnlySignal<TResult> Map<TSource, TResult>(this IReadOnlySignal<TSource> signal, Func<TSource, TResult> map)
        {
            return new CastingSignal<TSource, TResult>(signal, map);
        }

        struct CastingSignal<TSource, TResult> : IReadOnlySignal<TResult>
        {
            private IReadOnlySignal<TSource> _other;
            private Func<TSource, TResult> _map;

            public CastingSignal(IReadOnlySignal<TSource> other, Func<TSource, TResult> map)
            {
                _other = other ?? throw new ArgumentNullException(nameof(other));
                _map = map ?? throw new ArgumentNullException(nameof(map));
            }

            public IEnumerator<TResult> GetEnumerator()
            {
                return _other.Select(_map).GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return ((IEnumerable)_other).GetEnumerator();
            }

            public int Count => _other.Count;

            public TResult this[int index] => _map(_other[index]);

            public double Start => _other.Start;

            public double SamplingRate => _other.SamplingRate;

            public double SamplingPeriod => _other.SamplingPeriod;

            public double End => _other.End;
        }
        #endregion
    }
}
