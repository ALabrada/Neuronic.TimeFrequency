using System;
using System.Collections;
using System.Collections.Generic;

namespace Neuronic.TimeFrequency
{
    class ReadOnlySignal<T>: IReadOnlySignal<T>, IList<T>
    {
        public ReadOnlySignal(IList<T> samples, double start = 0d, double fs = 1d)
        {
            if (fs <= 0) throw new ArgumentOutOfRangeException(nameof(fs));
            Samples = samples ?? throw new ArgumentNullException(nameof(samples));
            Start = start;
            SamplingRate = fs;
            SamplingPeriod = 1d / fs;
            End = Start + (samples.Count == 0 ? 0 : SamplingPeriod * (samples.Count - 1));
        }

        public ReadOnlySignal(Signal<T> other)
        {
            Samples = other.Samples;
            Start = other.Start;
            SamplingRate = other.SamplingRate;
            SamplingPeriod = other.SamplingPeriod;
            End = other.End;
        }

        public double Start { get; }

        public double SamplingRate { get; }

        public double SamplingPeriod { get; }

        public double End { get; }

        public IList<T> Samples { get; }

        #region IEnumerable
        public IEnumerator<T> GetEnumerator()
        {
            return Samples.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return ((IEnumerable)Samples).GetEnumerator();
        }
        #endregion

        #region IList
        public void Add(T item)
        {
            throw new NotSupportedException();
        }

        public void Clear()
        {
            throw new NotSupportedException();
        }

        public bool Contains(T item)
        {
            return Samples.Contains(item);
        }

        public void CopyTo(T[] array, int arrayIndex)
        {
            Samples.CopyTo(array, arrayIndex);
        }

        public bool Remove(T item)
        {
            throw new NotSupportedException();
        }

        public int Count => Samples.Count;

        public bool IsReadOnly => true;

        public int IndexOf(T item)
        {
            return Samples.IndexOf(item);
        }

        public void Insert(int index, T item)
        {
            throw new NotSupportedException();
        }

        public void RemoveAt(int index)
        {
            throw new NotSupportedException();
        }

        public T this[int index]
        {
            get => Samples[index];
            set => throw new NotSupportedException();
        } 
        #endregion
    }
}