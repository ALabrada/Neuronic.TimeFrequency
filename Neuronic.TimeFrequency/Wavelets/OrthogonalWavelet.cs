using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.Interpolation;


namespace Neuronic.TimeFrequency.Wavelets
{
    /// <summary>
    /// Represents a wavelet defined by an orthogonal base.
    /// </summary>
    /// <seealso cref="Neuronic.TimeFrequency.Wavelets.WaveletBase" />
    public class OrthogonalWavelet : WaveletBase, IWavelet<double>
    {
        /// <summary>
        /// The precision level for <see cref="EvaluateDomain"/>.
        /// </summary>
        protected const int PrecisionLevel = 10;

        private readonly double[] _lowReconstructionFilter;
        private readonly double[] _highReconstructionFilter;
        private readonly double[] _lowDecompositionFilter;
        private readonly double[] _highDecompositionFilter;
        private readonly object _psiLock = new object();
        private Signal<double>? _psi;
        private double[] _x;

        /// <summary>
        /// Initializes a new instance of the <see cref="OrthogonalWavelet"/> class.
        /// </summary>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="lowRec">The low pass reconstruction filter.</param>
        /// <param name="highRec">The high pass reconstruction filter.</param>
        /// <param name="lowDec">The low pass decomposition filter.</param>
        /// <param name="highDec">The high pass decomposition filter.</param>
        /// <param name="vanishingMoments">The vanishing moments.</param>
        /// <param name="freq">The central frequency.</param>
        /// <exception cref="ArgumentNullException">Thrown when any of the filters is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">Thrown when the filter lengths do not match.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="vanishingMoments"/> is negative.</exception>
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
            if (vanishingMoments < 0) throw new ArgumentOutOfRangeException(nameof(vanishingMoments));

            base.CentralFrequency = freq;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="OrthogonalWavelet"/> class.
        /// </summary>
        /// <param name="shortName">The short name.</param>
        /// <param name="familyName">The family name.</param>
        /// <param name="filterLength">Length of the filter.</param>
        /// <param name="vanishingMoments">The vanishing moments.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="vanishingMoments"/> is negative.</exception>
        public OrthogonalWavelet(string shortName, string familyName, int filterLength, int vanishingMoments)
            : this (shortName, familyName, new double[filterLength], new double[filterLength], new double[filterLength], new double[filterLength], vanishingMoments)
        {
        }

        /// <summary>
        /// Gets the length of the filters.
        /// </summary>
        public int FilterLength { get; }

        /// <summary>
        /// Gets the amount of vanishing moments in the PSI function.
        /// </summary>
        public int VanishingMoments { get; }

        /// <summary>
        /// Gets the low pass reconstruction filter.
        /// </summary>
        public ReadOnlyCollection<double> LowReconstructionFilter => Array.AsReadOnly(_lowReconstructionFilter);

        /// <summary>
        /// Gets the high pass reconstruction filter.
        /// </summary>
        public ReadOnlyCollection<double> HighReconstructionFilter => Array.AsReadOnly(_highReconstructionFilter);

        /// <summary>
        /// Gets the low pass decomposition filter.
        /// </summary>
        public ReadOnlyCollection<double> LowDecompositionFilter => Array.AsReadOnly(_lowDecompositionFilter);

        /// <summary>
        /// Gets the high pass decomposition filter.
        /// </summary>
        public ReadOnlyCollection<double> HighDecompositionFilter => Array.AsReadOnly(_highDecompositionFilter);

        /// <summary>
        /// Gets the size that should have the signal returned by <see cref="EvaluateDomain"/>.
        /// </summary>
        protected virtual int DesiredOutputSize => (FilterLength - 1) * (1 << PrecisionLevel) + 1;

        /// <summary>
        /// Performs upsampling convolution with the filters to the specified level of precision.
        /// </summary>
        /// <param name="level">The level.</param>
        /// <returns>The convolved samples.</returns>
        protected virtual double[] Upcoef(int level)
        {
            var coeffs = new [] { Math.Pow(Math.Sqrt(2), level) };
            return Upcoef(coeffs, level);
        }

        /// <summary>
        /// Convolves the specified coefficient vector with the filters to the specified level of precision.
        /// </summary>
        /// <param name="coeffs">The coefficients to convolve.</param>
        /// <param name="level">The level.</param>
        /// <returns>The convolved coefficients.</returns>
        protected virtual double[] Upcoef(double[] coeffs, int level)
        {
            for (int i = 0; i < level; i++)
            {
                int recLen = 2 * coeffs.Length + FilterLength - 2;
                var rec = new double[recLen];
                if (i > 0)
                    UpsumplingConvolution(coeffs, _lowReconstructionFilter, rec);
                else
                    UpsumplingConvolution(coeffs, _highReconstructionFilter, rec);
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

        Signal<double> IWavelet<double>.Evaluate(double min, double max, int count)
        {
            if (count <= 0 || max < min)
                return new Signal<double>(new double[0]);

            var values = new double[count];
            var signal = new Signal<double>(values, min, (count - 1) / (max - min));
            Evaluate(signal);
            return signal;
        }

        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        public override void Evaluate(Signal<Complex> signal)
        {
            var phi = ProtectedEvaluate();
            Interpolate(phi, signal);
        }

        /// <summary>
        /// Evaluates the wavelet function in the range defined by the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        public virtual void Evaluate(Signal<double> signal)
        {
            var phi = ProtectedEvaluate();
            Interpolate(phi, signal);
        }

        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>
        /// The evaluated values.
        /// </returns>
        public override Signal<Complex> EvaluateDomain()
        {
            var psi = ProtectedEvaluate();
            var outputLength = Math.Max(psi.Count + 2, DesiredOutputSize);

            var values = new Complex[outputLength];
            var offset = 1;
            for (int i = 0; i < psi.Count; i++)
                values[i + offset] = psi[i];
           
            return new Signal<Complex>(values, psi.Start - (offset * psi.SamplingPeriod), psi.SamplingRate);
        }

        Signal<double> IWavelet<double>.EvaluateDomain()
        {
            var psi = ProtectedEvaluate();
            var outputLength = Math.Max(psi.Count + 2, DesiredOutputSize);

            var values = new double[outputLength];
            var offset = 1;
            Array.Copy(psi.Samples, 0, values, offset, psi.Count);

            return new Signal<double>(values, psi.Start - (offset * psi.SamplingPeriod), psi.SamplingRate);
        }

        /// <summary>
        /// Evaluates the wavelet function in all it's domain.
        /// </summary>
        /// <returns>
        /// The evaluated values without zero padding.
        /// </returns>
        protected virtual Signal<double> ProtectedEvaluate()
        {
            lock (_psiLock)
            {
                if (_psi.HasValue)
                    return _psi.Value;

                var p = 1 << PrecisionLevel;

                var keepLength = Enumerable.Range(0, PrecisionLevel).Aggregate(1, (total, _) => 2 * total + (FilterLength - 2));

                var psi = Upcoef(PrecisionLevel);
                var offset = 0;
                if (psi.Length > keepLength)
                {
                    offset = (psi.Length - keepLength) / 2;
                    Array.Clear(psi, 0, offset);
                    Array.Clear(psi, offset + keepLength, psi.Length - offset - keepLength);
                }

                _psi = new Signal<double>(psi, (double)(offset + 1) / p, p);
                return _psi.Value; 
            }
        }

        private double[] EvaluateTimeFor<T>(Signal<T> signal)
        {
            var x = new double[signal.Count];
            var min = signal.Start;
            var step = signal.SamplingPeriod;
            for (int i = 0; i < x.Length; i++)
                x[i] = min + i * step;
            return x;
        }

        /// <summary>
        /// Interpolates the specified function.
        /// </summary>
        /// <param name="input">The input function.</param>
        /// <param name="result">The result.</param>
        protected virtual void Interpolate(Signal<double> input, Signal<Complex> result)
        {
            var x = _x ?? (_x = EvaluateTimeFor(input));
            var y = input.Samples;
            var count = result.Count;
            var step = result.SamplingPeriod;
            var min = result.Start;
            var interpolator = StepInterpolation.InterpolateSorted(x, y);

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i] = t >= x[0] && t <= x[x.Length - 1] ? interpolator.Interpolate(t) : 0d;
            }
        }

        /// <summary>
        /// Interpolates the specified function.
        /// </summary>
        /// <param name="input">The input function.</param>
        /// <param name="result">The result.</param>
        protected virtual void Interpolate(Signal<double> input, Signal<double> result)
        {
            var x = _x ?? (_x = EvaluateTimeFor(input));
            var y = input.Samples;
            var count = result.Count;
            var step = result.SamplingPeriod;
            var min = result.Start;
            var interpolator = StepInterpolation.InterpolateSorted(x, y);

            for (int i = 0; i < count; i++)
            {
                var t = min + i * step;
                result[i] = t >= x[0] && t <= x[x.Length - 1] ? interpolator.Interpolate(t) : 0d;
            }
        }
    }
}
