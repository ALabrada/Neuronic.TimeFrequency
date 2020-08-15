using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neuronic.TimeFrequency
{
    public struct Signal<T>
        : IReadOnlySignal<T>
    {
        public Signal(T[] samples, double delay = 0d, double fs = 1d)
        {
            Samples = samples ?? throw new ArgumentNullException(nameof(samples));
            Delay = delay;
            SamplingRate = fs;
            SamplingPeriod = 1d / fs;
            Duration = SamplingPeriod * samples.Length;
        }

        public T this[int index]
        {
            get => Samples[index];
            set => Samples[index] = value;            
        }

        public T[] Samples { get; }

        public double Delay { get; }

        public double SamplingRate { get; }

        public double SamplingPeriod { get; }

        public double Duration { get; }

        public int Count => Samples.Length;

        public IEnumerator<T> GetEnumerator()
        {
            return Samples.AsEnumerable<T>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        } 
    }
}
