using System;
using System.Linq;
using System.Numerics;
using Accord.Math;

namespace Neuronic.TimeFrequency
{
    class OrthogonalWavelet : WaveletBase
    {
        private readonly double[] _lowReconstructionFilter;
        private readonly double[] _highReconstructionFilter;
        private readonly double[] _lowDecompositionFilter;
        private readonly double[] _highDecompositionFilter;
        private double _centralFrequency;

        public OrthogonalWavelet(double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, double freq = 0)
        {
            _lowReconstructionFilter = lowRec ?? throw new ArgumentNullException(nameof(lowRec));
            _highReconstructionFilter = highRec ?? throw new ArgumentNullException(nameof(highRec));
            _lowDecompositionFilter = lowDec ?? throw new ArgumentNullException(nameof(lowDec));
            _highDecompositionFilter = highDec ?? throw new ArgumentNullException(nameof(highDec));

            FilterLength = lowDec.Length;
            if (lowRec.Length != FilterLength || highDec.Length != FilterLength || highRec.Length != FilterLength)
                throw new ArgumentException("Filter length mismatch");

            _centralFrequency = freq;
        }

        public OrthogonalWavelet(int filterLength)
            : this (new double[filterLength], new double[filterLength], new double[filterLength], new double[filterLength])
        {
        }

        public int FilterLength { get; }

        public override Complex Energy { get; }

        public override double CentralFrequency => _centralFrequency > 0 ? _centralFrequency : (_centralFrequency = EstimateCentralFrequency(Upcoef(10), 1d / 1024d));

        public double[] LowReconstructionFilter => _lowReconstructionFilter;

        public double[] HighReconstructionFilter => _highReconstructionFilter;

        public double[] LowDecompositionFilter => _lowDecompositionFilter;

        public double[] HighDecompositionFilter => _highDecompositionFilter;

        private double[] Upcoef(int level, bool recA = false)
        {
            var coeffs = new double[] { Math.Pow(Math.Sqrt(2), level) };
            int recLen = 2 * coeffs.Length + FilterLength - 2;
            var rec = new double[recLen];
            for (int i = 0; i < level; i++)
            {
                if (recA || i > 0)
                    UpsumplingConvolution(coeffs, LowReconstructionFilter, rec);
                else
                    UpsumplingConvolution(coeffs, HighReconstructionFilter, rec);
            }
            return rec;
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
            const int level = 10;
            var p = 1 << level;

            var outputLength = (FilterLength - 1) * p + 1;
            var keepLength = Enumerable.Range(0, level).Aggregate(1, (total, _) => total = 2 * total + (FilterLength - 2));
            if (outputLength - keepLength - 2 < 0)
                outputLength = keepLength + 2;
            var rightExtentLength = outputLength - keepLength - 1;

            var psi = Upcoef(level, false);
            var offset = Math.Max(0, (psi.Length - keepLength) / 2);
            Array.Clear(psi, 0, offset);
            Array.Clear(psi, offset + keepLength, psi.Length - offset - keepLength);

            var x = new double[psi.Length];
            if (x.Length > keepLength)
            for (int i = 0; i < x.Length; i++)
                x[i] = (double) (i - offset + 1) / p;

            Interpolate(min, max, x, psi, values, start, count);
        }

        private void Interpolate(double min, double max, double[] x, double[] y, Complex[] result, int start, int count)
        {           
            var step = (max - min) / (count - 1);

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i + start] = Tools.Interpolate1D(t, x, y, 0, 0);
            }
        }
    }
}
