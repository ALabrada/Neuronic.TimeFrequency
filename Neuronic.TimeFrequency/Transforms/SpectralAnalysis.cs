using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.TimeFrequency.Transforms
{
    /// <summary>
    /// Represents the spectral analysis of a signal, by analyzing it's principal components.  
    /// </summary>
    public class SpectralAnalysis: IReadOnlyList<SpectralAnalysis.MonocomponentSignal>, ITimeFrequencyRepresentation
    {
        private readonly List<MonocomponentSignal> _components;

        /// <summary>
        /// Initializes a new instance of the <see cref="SpectralAnalysis"/> class.
        /// </summary>
        /// <param name="components">The signal components.</param>
        /// <param name="samplingPeriod">The sampling period.</param>
        /// <param name="sampleCount">The samples count.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="components"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="samplingPeriod"/> is not positive.</exception>
        public SpectralAnalysis(IEnumerable<MonocomponentSignal> components, double samplingPeriod, int sampleCount)
        {
            if (components == null) throw new ArgumentNullException(nameof(components));
            if (samplingPeriod <= 0) throw new ArgumentOutOfRangeException(nameof(samplingPeriod));
            SamplingPeriod = samplingPeriod;
            SampleCount = sampleCount;
            _components = components.ToList();
        }

        /// <summary>
        /// Gets the sampling period.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the amount of samples in the time domain.
        /// </summary>
        public int SampleCount { get; }

        /// <summary>
        /// Gets the number of components.
        /// </summary>
        public int Count => _components.Count;

        /// <summary>
        /// Gets the component at the specified index.
        /// </summary>
        /// <param name="index">The index.</param>
        /// <returns>The component at <paramref name="index"/>.</returns>
        public MonocomponentSignal this[int index] => _components[index];

        double ITimeFrequencyRepresentation.this[int offset, double frequency]
        {
            get
            {
                const double delta = 1d;
                var amplitude = 0d;

                foreach (var component in _components)
                {
                    var freqIndex = offset - component.FrequencyOffset;
                    if (freqIndex >= 0 && freqIndex < component.Frequency.Count &&
                        Math.Abs(component.Frequency[freqIndex] - frequency) < delta)
                        amplitude += component.Amplitude[offset];
                }

                return amplitude;
            }
        }

        /// <summary>
        /// Calculates the spectrogram of the signal by sampling the component signals in the frequency domain.
        /// </summary>
        /// <param name="spectralResolution">The desired spectral resolution.</param>
        /// <returns>The spectrogram.</returns>
        public IBilinearTimeFrequencyRepresentation GetSpectrogram(double spectralResolution = 0)
        {
            if (spectralResolution <= 0)
                spectralResolution = 1d / (SampleCount * SamplingPeriod);

            var maxFrequency = 0.5 / SamplingPeriod;
            var frequencies = new double[(int) Math.Ceiling(maxFrequency / spectralResolution)];
            for (int i = 0; i < frequencies.Length; i++)
                frequencies[i] = i * spectralResolution;

            var amplitudes = new double[SampleCount, frequencies.Length];
            for (int offset = 0; offset < SampleCount; offset++)
            {
                foreach (var component in _components)
                {
                    var freqIndex = offset - component.FrequencyOffset;
                    if (freqIndex >= 0 && freqIndex < component.Frequency.Count)
                    {
                        var frequency = component.Frequency[freqIndex];
                        var j = (int) Math.Round(frequency / spectralResolution);
                        if (j >= 0 && j < frequencies.Length)
                            amplitudes[offset, j] += component.Amplitude[offset];
                    }
                }
            }

            return new TimeFrequencyDistribution(amplitudes, frequencies, SamplingPeriod);
        }

        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>
        /// An enumerator that can be used to iterate through the collection.
        /// </returns>
        public IEnumerator<MonocomponentSignal> GetEnumerator()
        {
            return _components.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return ((IEnumerable)_components).GetEnumerator();
        }

        /// <summary>
        /// Represents a signal that has a well defined instantaneous frequency in its whole domain. 
        /// </summary>
        public struct MonocomponentSignal
        {
            /// <summary>
            /// Initializes a new instance of the <see cref="MonocomponentSignal"/> struct.
            /// </summary>
            /// <param name="amplitude">The amplitude of the signal.</param>
            /// <param name="frequency">The frequency of the signal.</param>
            /// <exception cref="ArgumentNullException">
            /// Throw when <paramref name="amplitude"/> or <paramref name="frequency"/> are <c>null</c>.
            /// </exception>
            public MonocomponentSignal(IReadOnlyList<double> amplitude, IReadOnlyList<double> frequency)
            {
                Amplitude = amplitude ?? throw new ArgumentNullException(nameof(amplitude));
                Frequency = frequency ?? throw new ArgumentNullException(nameof(frequency));
            }

            /// <summary>
            /// Gets the amplitude of the signal.
            /// </summary>
            public IReadOnlyList<double> Amplitude { get; }

            /// <summary>
            /// Gets the instantaneous frequency of the signal.
            /// </summary>
            public IReadOnlyList<double> Frequency { get; }

            /// <summary>
            /// Gets the frequency offset.
            /// </summary>
            /// <value>
            /// The frequency offset.
            /// </value>
            public int FrequencyOffset => (Amplitude.Count - Frequency.Count) / 2;
        }
    }
}