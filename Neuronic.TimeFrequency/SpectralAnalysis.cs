using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.TimeFrequency
{
    public class SpectralAnalysis: IReadOnlyList<SpectralAnalysis.MonocomponentSignal>
    {
        private readonly List<MonocomponentSignal> _components;

        public SpectralAnalysis(IEnumerable<MonocomponentSignal> components, double samplingPeriod)
        {
            if (components == null) throw new ArgumentNullException(nameof(components));
            if (samplingPeriod <= 0) throw new ArgumentOutOfRangeException(nameof(samplingPeriod));
            SamplingPeriod = samplingPeriod;
            _components = components.ToList();
        }

        public double SamplingPeriod { get; }

        public int Count => _components.Count;

        public MonocomponentSignal this[int index] => _components[index];

        public IEnumerator<MonocomponentSignal> GetEnumerator()
        {
            return _components.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return ((IEnumerable)_components).GetEnumerator();
        }

        public struct MonocomponentSignal
        {
            public MonocomponentSignal(IReadOnlyList<double> amplitude, IReadOnlyList<double> frequency)
            {
                Amplitude = amplitude ?? throw new ArgumentNullException(nameof(amplitude));
                Frequency = frequency ?? throw new ArgumentNullException(nameof(frequency));
            }

            public IReadOnlyList<double> Amplitude { get; }

            public IReadOnlyList<double> Frequency { get; }
        }
    }
}