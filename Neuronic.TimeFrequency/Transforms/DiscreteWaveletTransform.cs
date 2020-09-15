using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Represents the Discrete Wavelet Transform (DWT).  
    /// </summary>
    public class DiscreteWaveletTransform
    {
        private DiscreteWaveletTransform(IList<double> approximation, IList<double> detail, 
            double startTime, double samplingPeriod, OrthogonalWavelet wavelet)
        {
            Detail = new ReadOnlyCollection<double>(detail);
            Approximation = new ReadOnlyCollection<double>(approximation);
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
            StartTime = startTime;
        }

        private static DiscreteWaveletTransform Estimate(List<double> samples, double startTime, double samplingPeriod, OrthogonalWavelet wavelet)
        {
            var lowDec = wavelet.LowDecompositionFilter;
            var highDec = wavelet.HighDecompositionFilter;

            var length = (samples.Count - Math.Max(0, wavelet.FilterLength - 1)) / 2;
            var a = new double[length + 1];
            samples.Convolve(lowDec, a, 2);

            var d = new double[length + 1];
            samples.Convolve(highDec, d, 2);
            return new DiscreteWaveletTransform(
                new ArraySegment<double>(a, 1, length),
                new ArraySegment<double>(d, 1, length), startTime, samplingPeriod, wavelet);
        }

        /// <summary>
        /// Estimates the DWT for the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="padding">The padding method.</param>
        /// <returns>The computed DWT.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is <c>null</c>.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static DiscreteWaveletTransform Estimate(IReadOnlySignal<double> signal, OrthogonalWavelet wavelet, IPadding<double> padding)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            wavelet = wavelet ?? Wavelets.Wavelets.Haar;
            if (wavelet is BiorthogonalWavelet bior)
                wavelet = bior.Other;
            padding = padding ?? new SymetricPadding<double>();
            var samplingPeriod = signal.SamplingPeriod;

            var count = signal.Count + 2 * (wavelet.FilterLength - 1);
            var samples = new List<double>(count);
            samples.AddRange(padding.Pad(signal, wavelet.FilterLength - 1).Take(count));

            return Estimate(samples, signal.Start, samplingPeriod, wavelet);
        }

        /// <summary>
        /// Estimates the DWT for the specified signal.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="wavelet">The wavelet.</param>
        /// <param name="padding">The padding method.</param>
        /// <returns>The computed DWT.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is <c>null</c>.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>cwt</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        public static DiscreteWaveletTransform Estimate(IReadOnlySignal<float> signal, OrthogonalWavelet wavelet, IPadding<float> padding)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            wavelet = wavelet ?? Wavelets.Wavelets.Haar;
            if (wavelet is BiorthogonalWavelet bior)
                wavelet = bior.Other;
            padding = padding ?? new SymetricPadding<float>();
            var samplingPeriod = signal.SamplingPeriod;

            var count = signal.Count + 2 * (wavelet.FilterLength - 1);
            var samples = new List<double>(count);
            samples.AddRange(padding.Pad(signal, wavelet.FilterLength - 1).Take(count).Select(x => (double) x));

            return Estimate(samples, signal.Start, samplingPeriod, wavelet);
        }

        /// <summary>
        /// Gets the upper scale in a multi-scale DWT.
        /// </summary>
        public DiscreteWaveletTransform UpperScale { get; private set; }

        /// <summary>
        /// Gets the offset of the first sample in the time domain.
        /// </summary>
        double StartTime { get; }

        /// <summary>
        /// Gets the sampling period of the source signal.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the wavelet function.
        /// </summary>
        public OrthogonalWavelet Wavelet { get; }

        /// <summary>
        /// Gets the approximation coefficients.
        /// </summary>
        public ReadOnlyCollection<double> Approximation { get; }

        /// <summary>
        /// Gets the detail coefficients.
        /// </summary>
        public ReadOnlyCollection<double> Detail { get; }

        /// <summary>
        /// Estimates the next DWT scale.
        /// </summary>
        /// <param name="padding">The padding.</param>
        /// <returns>The next scale.</returns>
        public DiscreteWaveletTransform EstimateNextScale(IPadding<double> padding)
        {
            var signal = new ReadOnlySignal<double>(Approximation, StartTime, 2d / SamplingPeriod);
            var dwt = Estimate(signal, Wavelet, padding);
            dwt.UpperScale = this;
            return dwt;
        }

        /// <summary>
        /// Estimates a multi-scale DWT analysis of the signal.
        /// </summary>
        /// <param name="padding">The padding.</param>
        /// <param name="scales">The maximum number of scales or <c>zero</c> for unlimited.</param>
        /// <returns>The lowest scale of the multi-scale DWT.</returns>
        public DiscreteWaveletTransform EstimateMultiscale(IPadding<double> padding, int scales = -1)
        {
            if (Approximation.Count < 2 || scales == 0)
                return this;
            var next = EstimateNextScale(padding);
            return next.EstimateMultiscale(padding, scales - 1);
        }

        /// <summary>
        /// Reconstructs the original signal by reversing the DWT.
        /// </summary>
        /// <returns>The reconstructed signal.</returns>
        /// <remarks>
        /// If the DWT is multi-scale, the reconstruction includes also the upper scales.
        /// </remarks>
        public IReadOnlySignal<double> Reverse()
        {
            var samples = ReverseMultiscale();
            return new ReadOnlySignal<double>(samples, StartTime, 1d / SamplingPeriod);
        }

        /// <summary>
        /// Reverses the DWT using the approximation coefficients reconstructed from lower scales.
        /// </summary>
        /// <param name="approximation">The approximation coefficients.</param>
        /// <returns>The reconstructed signal.</returns>
        protected IList<double> ReverseMultiscale(IList<double> approximation = null)
        {
            approximation = approximation ?? Approximation;
            var signal = new double[Detail.Count * 2 + Wavelet.FilterLength - 2];
            var lowRec = Wavelet.LowReconstructionFilter;
            var highRec = Wavelet.HighReconstructionFilter;

            OrthogonalWavelet.UpsumplingConvolution(approximation, lowRec, signal);
            OrthogonalWavelet.UpsumplingConvolution(Detail, highRec, signal);

            var result = new ArraySegment<double>(signal, Wavelet.FilterLength - 2, 2 * Detail.Count - Wavelet.FilterLength + 2) as IList<double>;
            if (UpperScale != null)
                result = UpperScale.ReverseMultiscale(result);

            return result;
        }
    }
}