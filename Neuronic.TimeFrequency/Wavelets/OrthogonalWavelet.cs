using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Accord.Math;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class OrthogonalWavelet : WaveletBase
    {
        protected const int PrecissionLevel = 10;

        private readonly double[] _lowReconstructionFilter;
        private readonly double[] _highReconstructionFilter;
        private readonly double[] _lowDecompositionFilter;
        private readonly double[] _highDecompositionFilter;
        private readonly Dictionary<int, double[]> _psiCache = new Dictionary<int, double[]>();

        public OrthogonalWavelet(string shortName, string familyName, 
            double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, 
            int vanishingMoments, double freq = 0)
        : base (shortName, familyName)
        {
            VanishingMoments = vanishingMoments;
            _lowReconstructionFilter = lowRec ?? throw new ArgumentNullException(nameof(lowRec));
            _highReconstructionFilter = highRec ?? throw new ArgumentNullException(nameof(highRec));
            _lowDecompositionFilter = lowDec ?? throw new ArgumentNullException(nameof(lowDec));
            _highDecompositionFilter = highDec ?? throw new ArgumentNullException(nameof(highDec));

            FilterLength = lowDec.Length;
            if (lowRec.Length != FilterLength || highDec.Length != FilterLength || highRec.Length != FilterLength)
                throw new ArgumentException("Filter length mismatch");

            base.CentralFrequency = freq;
        }

        public OrthogonalWavelet(string shortName, string familyName, int filterLength, int vanishingMoments)
            : this (shortName, familyName, new double[filterLength], new double[filterLength], new double[filterLength], new double[filterLength], vanishingMoments)
        {
        }

        public int FilterLength { get; }

        public int VanishingMoments { get; }

        public double[] LowReconstructionFilter => _lowReconstructionFilter;

        public double[] HighReconstructionFilter => _highReconstructionFilter;

        public double[] LowDecompositionFilter => _lowDecompositionFilter;

        public double[] HighDecompositionFilter => _highDecompositionFilter;

        protected virtual int DesiredOutputSize => (FilterLength - 1) * (1 << PrecissionLevel) + 1;

        protected virtual double[] Upcoef(int level)
        {
            if (_psiCache.TryGetValue(level, out var psi))
                return psi;
            var coeffs = new [] { Math.Pow(Math.Sqrt(2), level) };
            psi = Upcoef(coeffs, level);
            _psiCache[level] = psi;
            return psi;
        }

        protected virtual double[] Upcoef(double[] coeffs, int level)
        {
            for (int i = 0; i < level; i++)
            {
                int recLen = 2 * coeffs.Length + FilterLength - 2;
                var rec = new double[recLen];
                if (i > 0)
                    UpsumplingConvolution(coeffs, LowReconstructionFilter, rec);
                else
                    UpsumplingConvolution(coeffs, HighReconstructionFilter, rec);
                coeffs = rec;
            }

            return coeffs;
        }

        private void UpsumplingConvolution(double[] coeff, double[] filter, double[] output)
        {
            if (filter.Length < 2 || filter.Length % 2 != 0)
                throw new ArgumentException("Invalid filter length", nameof(filter));

            var i = 0;
            var o = 0;

            for (; i < coeff.Length && i < filter.Length / 2; ++i, o += 2)
            {
                for (var j = 0; j <= i; ++j)
                {
                    output[o] += filter[j * 2] * coeff[i - j];
                    output[o + 1] += filter[j * 2 + 1] * coeff[i - j];
                }
            }

            for (; i < coeff.Length; ++i, o += 2)
            {
                for (var j = 0; j < filter.Length / 2; ++j)
                {
                    output[o] += filter[j * 2] * coeff[i - j];
                    output[o + 1] += filter[j * 2 + 1] * coeff[i - j];
                }
            }

            for (; i < filter.Length / 2; ++i, o += 2)
            {
                for (var j = i - (coeff.Length - 1); j <= i; ++j)
                {
                    output[o] += filter[j * 2] * coeff[i - j];
                    output[o + 1] += filter[j * 2 + 1] * coeff[i - j];
                }
            }

            for (; i < coeff.Length + filter.Length / 2; ++i, o += 2)
            {
                for (var j = i - (coeff.Length - 1); j < filter.Length / 2; ++j)
                {
                    output[o] += filter[j * 2] * coeff[i - j];
                    output[o + 1] += filter[j * 2 + 1] * coeff[i - j];
                }
            }
        }

        public override void Evaluate(Signal<Complex> signal)
        {
            var phi = ProtectedEvaluate();
            Interpolate(phi, signal);
        }

        public virtual void Evaluate(Signal<double> signal)
        {
            var phi = ProtectedEvaluate();
            Interpolate(phi, signal);
        }

        public override Signal<Complex> Evaluate()
        {
            var psi = ProtectedEvaluate();
            var outputLength = Math.Max(psi.Count + 2, DesiredOutputSize);

            var values = new Complex[outputLength];
            var offset = 1;
            for (int i = 0; i < psi.Count; i++)
                values[i + offset] = psi[i];
           
            return new Signal<Complex>(values, psi.Delay - (offset * psi.SamplingPeriod), psi.SamplingRate);
        }

        protected virtual Signal<double> ProtectedEvaluate()
        {
            var p = 1 << PrecissionLevel;

            var keepLength = Enumerable.Range(0, PrecissionLevel).Aggregate(1, (total, _) => 2 * total + (FilterLength - 2));

            var psi = Upcoef(PrecissionLevel);
            var offset = 0;
            if (psi.Length > keepLength)
            {
                offset = (psi.Length - keepLength) / 2;
                Array.Clear(psi, 0, offset);
                Array.Clear(psi, offset + keepLength, psi.Length - offset - keepLength);
            }

            return new Signal<double>(psi, (double)(offset + 1) / p, p);
        }

        protected virtual double[] EvaluateTimeFor<T>(Signal<T> signal)
        {
            var x = new double[signal.Count];
            var min = signal.Delay;
            var step = signal.SamplingPeriod;
            for (int i = 0; i < x.Length; i++)
                x[i] = min + i * step;
            return x;
        }

        protected virtual void Interpolate(Signal<double> input, Signal<Complex> result)
        {
            var x = EvaluateTimeFor(input);
            var y = input.Samples;
            var count = result.Count;
            var step = result.SamplingPeriod;
            var min = result.Delay;

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i] = Tools.Interpolate1D(t, x, y, 0, 0);
                //var index = Array.BinarySearch(x, t);
                //if (index >= 0)
                //    result[i + start] = y[index];
                //else
                //{
                //    index = ~index;
                //    if (index == 0 || index == x.Length)
                //        result[i] = 0d;
                //    else
                //        result[i] = y[index - 1];
                //}
            }
        }

        protected virtual void Interpolate(Signal<double> input, Signal<double> result)
        {
            var x = EvaluateTimeFor(input);
            var y = input.Samples;
            var count = result.Count;
            var step = result.SamplingPeriod;
            var min = result.Delay;

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i] = Tools.Interpolate1D(t, x, y, 0, 0);
                //var index = Array.BinarySearch(x, t);
                //if (index >= 0)
                //    result[i + start] = y[index];
                //else
                //{
                //    index = ~index;
                //    if (index == 0 || index == x.Length)
                //        result[i] = 0d;
                //    else
                //        result[i] = y[index - 1];
                //}
            }
        }
    }
}
