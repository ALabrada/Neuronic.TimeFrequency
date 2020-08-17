using System;
using System.Collections.Generic;

namespace Neuronic.TimeFrequency
{
    /// <summary>
    /// Abstraction of a signal padding method.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public interface IPadding<T>
    {
        /// <summary>
        /// Pads the specified list in both directions with the specified count of elements.
        /// </summary>
        /// <param name="list">The list.</param>
        /// <param name="count">The count.</param>
        /// <returns>The padded sequence.</returns>
        IEnumerable<T> Pad(IReadOnlyList<T> list, int count);
    }

    /// <summary>
    /// A padding method that reflects the list elements across it's boundaries.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <seealso cref="Neuronic.TimeFrequency.IPadding{T}" />
    public class SymetricPadding<T> : IPadding<T>
    {
        /// <inheritdoc />
        public IEnumerable<T> Pad(IReadOnlyList<T> list, int count)
        {
            if (count > list.Count)
                throw new ArgumentOutOfRangeException(nameof(count));
            for (int i = count - 1; i >= 0; i--)
                yield return list[i];
            foreach (var item in list)
                yield return item;
            for (int i = list.Count - 1; i >= list.Count - count; i--)
                yield return list[i];
        }
    }

    /// <summary>
    /// Simple padding method that uses zeros as padding elements.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <seealso cref="Neuronic.TimeFrequency.IPadding{T}" />
    public class ZeroPadding<T> : IPadding<T>
    {
        /// <inheritdoc />
        public IEnumerable<T> Pad(IReadOnlyList<T> list, int count)
        {
            for (int i = count - 1; i >= 0; i--)
                yield return default(T);
            foreach (var item in list)
                yield return item;
            for (int i = count - 1; i >= 0; i--)
                yield return default(T);
        }
    }
}