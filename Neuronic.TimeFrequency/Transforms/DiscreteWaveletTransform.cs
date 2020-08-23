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
        private DiscreteWaveletTransform(IList<double> approximation, IList<double> detail, double samplingPeriod, OrthogonalWavelet wavelet)
        {
            Detail = new ReadOnlyCollection<double>(detail);
            Approximation = new ReadOnlyCollection<double>(approximation);
            SamplingPeriod = samplingPeriod;
            Wavelet = wavelet;
        }

        private static DiscreteWaveletTransform Estimate(List<double> samples, double samplingPeriod, OrthogonalWavelet wavelet)
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
                new ArraySegment<double>(d, 1, length),
                samplingPeriod, wavelet);
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

            return Estimate(samples, samplingPeriod, wavelet);
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

            return Estimate(samples, samplingPeriod, wavelet);
        }

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
    }
}