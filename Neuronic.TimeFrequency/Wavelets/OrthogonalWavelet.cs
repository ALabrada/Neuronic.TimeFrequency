using System;
using System.Linq;
using System.Numerics;
using Accord.Math;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class OrthogonalWavelet : WaveletBase
    {
        private readonly double[] _lowReconstructionFilter;
        private readonly double[] _highReconstructionFilter;
        private readonly double[] _lowDecompositionFilter;
        private readonly double[] _highDecompositionFilter;
        private double _centralFrequency;

        public OrthogonalWavelet(double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, double freq = 0)
        {
            VanishingMoments = vanishingMoments;
            _lowReconstructionFilter = lowRec ?? throw new ArgumentNullException(nameof(lowRec));
            _highReconstructionFilter = highRec ?? throw new ArgumentNullException(nameof(highRec));
            _lowDecompositionFilter = lowDec ?? throw new ArgumentNullException(nameof(lowDec));
            _highDecompositionFilter = highDec ?? throw new ArgumentNullException(nameof(highDec));

            FilterLength = lowDec.Length;
            if (lowRec.Length != FilterLength || highDec.Length != FilterLength || highRec.Length != FilterLength)
                throw new ArgumentException("Filter length mismatch");

            _centralFrequency = freq;
        }

        public OrthogonalWavelet(int filterLength, int vanishingMoments)
            : this (new double[filterLength], new double[filterLength], new double[filterLength], new double[filterLength], vanishingMoments)
        {
        }

        public int FilterLength { get; }

        public int VanishingMoments { get; }

        public override Complex Energy { get; }

        public override double CentralFrequency => _centralFrequency > 0 ? _centralFrequency : (_centralFrequency = EstimateCentralFrequency(10));

        protected double EstimateCentralFrequency(int level)
        {
            var p = 1 << level;
            var phi = Upcoef(level);
            var values = new Complex[Tools.NextPowerOf2(phi.Length + 2)];
            var offset = (values.Length - phi.Length) / 2;
            for (int i = 0; i < phi.Length; i++)
                values[i + offset] = phi[i];
            return EstimateCentralFrequency(values, 1d / p);
        }

        public double[] LowReconstructionFilter => _lowReconstructionFilter;

        public double[] HighReconstructionFilter => _highReconstructionFilter;

        public double[] LowDecompositionFilter => _lowDecompositionFilter;

        public double[] HighDecompositionFilter => _highDecompositionFilter;

        protected virtual double[] Upcoef(int level, bool recA = false)
        {
            var coeffs = new [] { Math.Pow(Math.Sqrt(2), level) };
            return Upcoef(coeffs, level, recA);
        }

        protected virtual double[] Upcoef(double[] coeffs, int level, bool recA = false)
        {
            for (int i = 0; i < level; i++)
            {
                int recLen = 2 * coeffs.Length + FilterLength - 2;
                var rec = new double[recLen];
                if (recA || i > 0)
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

        public override void Evaluate(double min, double max, Complex[] values, int start, int count)
        {
            Evaluate(out var x, out var psi);
            Interpolate(min, max, x, psi, values, start, count);
        }

        protected virtual void Evaluate(out double[] x, out double[] psi)
        {
            const int level = 10;
            var p = 1 << level;

            var keepLength = Enumerable.Range(0, level).Aggregate(1, (total, _) => 2 * total + (FilterLength - 2));

            psi = Upcoef(level, false);
            var offset = 0;
            if (psi.Length > keepLength)
            {
                offset = (psi.Length - keepLength) / 2;
                Array.Clear(psi, 0, offset);
                Array.Clear(psi, offset + keepLength, psi.Length - offset - keepLength);
            }

            x = new double[psi.Length];
            for (int i = 0; i < x.Length; i++)
                x[i] = (double) (i - offset + 1) / p;
        }

        protected virtual void Interpolate(double min, double max, double[] x, double[] y, Complex[] result, int start, int count)
        {           
            var step = (max - min) / (count - 1);

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i + start] = Tools.Interpolate1D(t, x, y, 0, 0);
                //var index = Array.BinarySearch(x, t);
                //if (index >= 0)
                //    result[i + start] = y[index];
                //else
                //{
                //    index = ~index;
                //    if (index == 0 || index == x.Length)
                //        result[i + start] = 0d;
                //    else
                //        result[i + start] = y[index - 1];
                //}
            }
        }
    }
}
