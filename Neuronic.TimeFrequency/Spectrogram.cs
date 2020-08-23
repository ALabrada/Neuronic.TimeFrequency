using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Accord.Math;
using Accord.Math.Transforms;
using Neuronic.TimeFrequency.Kernels;

namespace Neuronic.TimeFrequency
{
    public class Spectrogram: IBilinearTimeFrequencyRepresentation
    {
        private Complex[,] _values;
        private double[] _frequencies;

        internal Spectrogram(Complex[,] values, double[] frequencies, double startTime, double samplingPeriod)
        {
            _values = values ?? throw new ArgumentNullException(nameof(values));
            _frequencies = frequencies ?? throw new ArgumentNullException(nameof(frequencies));
            SamplingPeriod = samplingPeriod;
            StartTime = startTime;
        }

        public static Spectrogram Estimate(IReadOnlySignal<double> signal, Func<int, double[]> winFunc = null, int? overlap = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            var window = winFunc?.Invoke(signal.Count) ?? Windows.Hamming.Evaluate(signal.Count);
            var stride = Math.Max(1, window.Length - (overlap ?? window.Length / 2));
            var sampCount = (signal.Count + stride - window.Length) / stride;

            var wFreq = 2 * Math.PI / signal.SamplingPeriod;
            var nfft = Math.Max(256, Tools.NextPowerOf2(window.Length));
            var freqCount = nfft / 2 + 1;
            var frequencies = new double[freqCount];
            for (int i = 0; i < freqCount; i++)
                frequencies[i] = wFreq * i / nfft;

            var values = new Complex[sampCount, freqCount];
            var windowed = new Complex[nfft];

            for (int offset = 0; offset < sampCount; offset++)
            {
                for (int i = 0, k = offset * stride; i < window.Length; i++, k++)
                    if (k < signal.Count)
                        windowed[i] = window[i] * signal[k];

                FourierTransform2.FFT(windowed, FourierTransform.Direction.Forward);

                for (int i = 0; i < freqCount; i++)
                    values[offset, i] = windowed[i];

                Array.Clear(windowed, 0, windowed.Length);
            }

            return new Spectrogram(values, frequencies, signal.Start + (window.Length / 2) * signal.SamplingPeriod, signal.SamplingPeriod * stride);
        }

        public static Spectrogram Estimate(IReadOnlySignal<double> signal, double[] window, int? overlap = null)
        {
            return Estimate(signal, _ => window, overlap);
        }

        public static Spectrogram Estimate(IReadOnlySignal<float> signal, Func<int, double[]> winFunc = null, int? overlap = null)
        {
            if (signal == null) throw new ArgumentNullException(nameof(signal));
            return Estimate(signal.Map(x => (double) x), winFunc, overlap);
        }

        public static Spectrogram Estimate(IReadOnlySignal<float> signal, double[] window, int? overlap = null)
        {
            return Estimate(signal, _ => window, overlap);
        }

        public double StartTime { get; }

        public double SamplingPeriod { get; }

        public int SampleCount => _values.GetLength(0);

        public int FrequencyCount => _values.GetLength(1);

        public IEnumerable<double> Frequencies => _frequencies;

        double ITimeFrequencyRepresentation.this[int offset, double frequency] =>
            this[offset, frequency].SquaredMagnitude();

        double IBilinearTimeFrequencyRepresentation.this[int offset, int frequencyIndex] =>
            this[offset, frequencyIndex].SquaredMagnitude();

        public Complex this[int offset, int freqIndex] => _values[offset, freqIndex];

        public Complex this[double delay, double frequency]
        {
            get
            {
                var offset = (int)Math.Round(delay / SamplingPeriod);
                var freqIndex = Array.BinarySearch(_frequencies, frequency);
                if (freqIndex < 0)
                    freqIndex = ~freqIndex;
                return this[offset, freqIndex];
            }
        }

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
        public IEnumerable<Complex> EnumerateFrequency(int index) =>
            Enumerable.Range(0, _values.GetLength(0)).Select(j => _values[j, index]);
    }
}
