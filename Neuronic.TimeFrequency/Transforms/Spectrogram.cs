using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using MathNet.Numerics;
using Neuronic.TimeFrequency.Kernels;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Spectrogram using Short Time Fourier Transform (STFT).
    /// </summary>
    public class Spectrogram: IBilinearTimeFrequencyRepresentation
    {
        private readonly Complex[,] _values;
        private readonly double[] _frequencies;

        internal Spectrogram(Complex[,] values, double[] frequencies, double startTime, double samplingPeriod)
        {
            _values = values ?? throw new ArgumentNullException(nameof(values));
            _frequencies = frequencies ?? throw new ArgumentNullException(nameof(frequencies));
            SamplingPeriod = samplingPeriod;
            StartTime = startTime;
        }

        /// <summary>
        /// Estimates the Spectrogram of the specified signal using the STFT algorithm.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="winFunc">The window function. Default to hamming window.</param>
        /// <param name="overlap">The number of overlapped samples between consecutive window offsets.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The Spectrogram of <paramref name="signal"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is null.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>spectrogram</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        /// <seealso cref="Window"/>
        public static Spectrogram Estimate(IReadOnlySignal<double> signal, Func<int, double[]> winFunc = null, int? overlap = null, 
            ParallelOptions options = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            options = options ?? new ParallelOptions();
            var window = winFunc?.Invoke(signal.Count) ?? Window.Hamming(signal.Count);
            var stride = Math.Max(1, window.Length - (overlap ?? window.Length / 2));
            var sampCount = (signal.Count + stride - window.Length) / stride;

            var wFreq = 2 * Math.PI / signal.SamplingPeriod;
            var nfft = Math.Max(256, window.Length.NextPowerOf2());
            var freqCount = nfft / 2 + 1;
            var frequencies = new double[freqCount];
            for (int i = 0; i < freqCount; i++)
                frequencies[i] = wFreq * i / nfft;

            var values = new Complex[sampCount, freqCount];

            Parallel.For(0, sampCount, options,
                () => new { windowed = new Complex[nfft] },
                (offset, _, state) =>
                {
                    var windowed = state.windowed;

                    for (int i = 0, k = offset * stride; i < window.Length; i++, k++)
                        if (k < signal.Count)
                            windowed[i] = window[i] * signal[k];

                    windowed.FFT();

                    for (int i = 0; i < freqCount; i++)
                        values[offset, i] = windowed[i];

                    Array.Clear(windowed, 0, windowed.Length);

                    return state;
                },
                _ => {});

            return new Spectrogram(values, frequencies, signal.Start + (window.Length / 2) * signal.SamplingPeriod, signal.SamplingPeriod * stride);
        }

        /// <summary>
        /// Estimates the Spectrogram of the specified signal using the STFT algorithm.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="window">The window. Default to hamming window.</param>
        /// <param name="overlap">The number of overlapped samples between consecutive window offsets.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The Spectrogram of <paramref name="signal"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is null.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>spectrogram</c> function in <c>Matlab R2014</c>.
        /// </remarks> 
        /// <seealso cref="Window"/>
        public static Spectrogram Estimate(IReadOnlySignal<double> signal, double[] window, int? overlap = null,
            ParallelOptions options = null)
        {
            return Estimate(signal, _ => window, overlap, options);
        }

        /// <summary>
        /// Estimates the Spectrogram of the specified signal using the STFT algorithm.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="winFunc">The window function. Default to hamming window.</param>
        /// <param name="overlap">The number of overlapped samples between consecutive window offsets.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The Spectrogram of <paramref name="signal"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is null.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>spectrogram</c> function in <c>Matlab R2014</c>.
        /// </remarks> 
        /// <seealso cref="Window"/>
        public static Spectrogram Estimate(IReadOnlySignal<float> signal, Func<int, double[]> winFunc = null, int? overlap = null,
            ParallelOptions options = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            return Estimate(signal.Map(x => (double) x), winFunc, overlap, options);
        }

        /// <summary>
        /// Estimates the Spectrogram of the specified signal using the STFT algorithm.
        /// </summary>
        /// <param name="signal">The signal.</param>
        /// <param name="window">The window. Default to hamming window.</param>
        /// <param name="overlap">The number of overlapped samples between consecutive window offsets.</param>
        /// <param name="options">The options for parallelization.</param>
        /// <returns>The Spectrogram of <paramref name="signal"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="signal"/> is null.</exception>
        /// <remarks>
        /// This algorithm is based on the <c>spectrogram</c> function in <c>Matlab R2014</c>.
        /// </remarks>
        /// <seealso cref="Window"/>
        public static Spectrogram Estimate(IReadOnlySignal<float> signal, double[] window, int? overlap = null,
            ParallelOptions options = null)
        {
            return Estimate(signal, _ => window, overlap, options);
        }

        /// <summary>
        /// Gets the offset of the first sample in the time domain.
        /// </summary>
        public double StartTime { get; }

        /// <summary>
        /// Gets the sampling period in the time domain.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of samples in the time domain.
        /// </summary>
        public int SampleCount => _values.GetLength(0);

        /// <summary>
        /// Gets the amount of samples in the frequency domain.
        /// </summary>
        public int FrequencyCount => _values.GetLength(1);

        /// <summary>
        /// Gets the computed frequencies.
        /// </summary>
        public IEnumerable<double> Frequencies => _frequencies;

        double ITimeFrequencyRepresentation.this[int offset, double frequency] =>
            this[offset, frequency].MagnitudeSquared();

        double IBilinearTimeFrequencyRepresentation.this[int offset, int frequencyIndex] =>
            this[offset, frequencyIndex].MagnitudeSquared();

        /// <summary>
        /// Gets the spectrogram value for the specified offset and frequency.
        /// </summary>
        /// <param name="offset">The sample offset in time domain.</param>
        /// <param name="freqIndex">The sample offset in the frequency domain.</param>
        /// <returns>The time-frequency content.</returns>
        public Complex this[int offset, int freqIndex] => _values[offset, freqIndex];

        /// <summary>
        /// Gets the spectrogram value for the specified offset and frequency.
        /// </summary>
        /// <param name="time">The time.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The time-frequency content.</returns>
        public Complex this[double time, double frequency]
        {
            get
            {
                var offset = this.FindClosestOffsetOfTime(time);
                var freqIndex = this.FindClosestIndexOfFrequency(frequency);
                return this[offset, freqIndex];
            }
        }

        /// <summary>
        /// Gets the spectrogram value for the specified offset and frequency.
        /// </summary>
        /// <param name="offset">The sample offset in time domain.</param>
        /// <param name="frequency">The frequency.</param>
        /// <returns>The time-frequency content.</returns>
        public Complex this[int offset, double frequency]
        {
            get
            {
                var freqIndex = Array.BinarySearch(_frequencies, frequency);
                if (freqIndex < 0)
                    freqIndex = ~freqIndex;
                return this[offset, freqIndex];
            }
        }

        /// <summary>
        /// Enumerates the values associated with the specified frequency.
        /// </summary>
        /// <param name="index">The index of the frequency.</param>
        /// <returns>The TFD values.</returns>
        public IEnumerable<Complex> EnumerateValuesOfFrequencyAt(int index) =>
            Enumerable.Range(0, _values.GetLength(0)).Select(j => _values[j, index]);
    }
}
