using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Neuronic.TimeFrequency
{
    /// <summary>
    /// Represents a measured signal sampled with a fixed rate.
    /// </summary>
    /// <typeparam name="T">The type of the samples.</typeparam>
    /// <seealso cref="Neuronic.TimeFrequency.IReadOnlySignal{T}" />
    public struct Signal<T> : IList<T>, 
            IReadOnlySignal<T>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Signal{T}"/> type.
        /// </summary>
        /// <param name="samples">The samples.</param>
        /// <param name="start">The time of the first sample.</param>
        /// <param name="fs">The sampling rate.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="samples"/> is <c>null</c>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if <paramref name="fs"/> is not a positive value.</exception>
        public Signal(T[] samples, double start = 0d, double fs = 1d)
        {
            if (fs <= 0) throw new ArgumentOutOfRangeException(nameof(fs));
            Samples = samples ?? throw new ArgumentNullException(nameof(samples));
            Start = start;
            SamplingRate = fs;
            SamplingPeriod = 1d / fs;
            End = Start + (samples.Length == 0 ? 0 :SamplingPeriod * (samples.Length - 1));
        }

        /// <summary>
        /// Gets or sets the sample at the specified index.
        /// </summary>
        /// <param name="index">The index.</param>
        /// <returns>The sample at <paramref name="index"/>.</returns>
        public T this[int index]
        {
            get => Samples[index];
            set => Samples[index] = value;            
        }

        /// <summary>
        /// Gets the samples.
        /// </summary>
        internal T[] Samples { get; }

        /// <summary>
        /// Gets the time of the first sample.
        /// </summary>
        public double Start { get; }

        /// <summary>
        /// Gets the sampling rate.
        /// </summary>
        public double SamplingRate { get; }

        /// <summary>
        /// Gets the sampling period.
        /// </summary>
        public double SamplingPeriod { get; }

        /// <summary>
        /// Gets the time of the last sample.
        /// </summary>
        public double End { get; }

        /// <summary>
        /// Gets the amount of samples.
        /// </summary>
        public int Count => Samples.Length;

        #region IEnumerable
        /// <summary>
        /// Returns an enumerator that iterates through the collection.
        /// </summary>
        /// <returns>
        /// An enumerator that can be used to iterate through the collection.
        /// </returns>
        public IEnumerator<T> GetEnumerator()
        {
            return Samples.AsEnumerable<T>().GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        } 
        #endregion

        #region IList
        int IList<T>.IndexOf(T item)
        {
            return Array.IndexOf(Samples, item);
        }

        void IList<T>.Insert(int index, T item)
        {
            throw new NotSupportedException();
        }

        void IList<T>.RemoveAt(int index)
        {
            throw new NotSupportedException();
        }

        void ICollection<T>.Add(T item)
        {
            throw new NotSupportedException();
        }

        void ICollection<T>.Clear()
        {
            throw new NotSupportedException();
        }

        bool ICollection<T>.Contains(T item)
        {
            return Array.IndexOf(Samples, item) >= 0;
        }

        void ICollection<T>.CopyTo(T[] array, int arrayIndex)
        {
            Array.Copy(Samples, 0, array, arrayIndex, Count);
        }

        bool ICollection<T>.Remove(T item)
        {
            throw new NotSupportedException();
        }

        bool ICollection<T>.IsReadOnly => Samples.IsReadOnly;
        #endregion
    }
}
