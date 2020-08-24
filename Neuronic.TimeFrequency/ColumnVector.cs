using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Neuronic.TimeFrequency
{
    class ColumnVector<T>: IReadOnlyList<T>
    {
        private readonly T[,] _matrix;

        public ColumnVector(T[,] matrix, int columnIndex)
        {
            _matrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (columnIndex < 0 || columnIndex >= _matrix.GetLength(1)) throw new ArgumentOutOfRangeException(nameof(columnIndex));
            ColumnIndex = columnIndex;
        }

        public int ColumnIndex { get; }

        public int Count => _matrix.GetLength(0);

        public T this[int index] => _matrix[index, ColumnIndex];

        public IEnumerator<T> GetEnumerator()
        {
            return Enumerable.Range(0, Count).Select(i => this[i]).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    class RowVector<T> : IReadOnlyList<T>
    {
        private readonly T[,] _matrix;

        public RowVector(T[,] matrix, int rowIndex)
        {
            _matrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (rowIndex < 0 || rowIndex >= _matrix.GetLength(0)) throw new ArgumentOutOfRangeException(nameof(rowIndex));
            RowIndex = rowIndex;
        }

        public int RowIndex { get; }

        public int Count => _matrix.GetLength(1);

        public T this[int index] => _matrix[RowIndex, index];

        public IEnumerator<T> GetEnumerator()
        {
            return Enumerable.Range(0, Count).Select(i => this[i]).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}