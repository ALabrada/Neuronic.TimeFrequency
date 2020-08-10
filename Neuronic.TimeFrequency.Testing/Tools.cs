using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Neuronic.TimeFrequency.Testing
{
    static class Tools
    {
        public static IList<float> ReadColumn(this TextReader reader)
        {
            var values = new List<float>();
            string line;
            while ((line = reader.ReadLine()) != null)
                values.AddRange(ReadNumbersFrom(line));
            return values;
        }

        private static char[] Separators = new [] { ' ', ',', '\n', '\r', '\t' };
        public static IEnumerable<float> ReadNumbersFrom(string str)
        {
            var parts = str.Split(Separators, StringSplitOptions.RemoveEmptyEntries);
            foreach (var part in parts)
            {
                if (float.TryParse(part, NumberStyles.Any, CultureInfo.InvariantCulture, out var value))
                    yield return value;
            }
        }

        public static void AssertAreEqual(IList<Complex> expectedValues, IList<Complex> actualValues, double delta)
        {
            Assert.AreEqual(expectedValues.Count, actualValues.Count);
            var dif = expectedValues.Zip(actualValues, (x, y) => (x - y).Magnitude).Sum();
            var relDif = dif / expectedValues.Count;
            Assert.AreEqual(0, relDif, delta);
        }

        public static void AssertAreEqual(IList<double> expectedValues, IList<double> actualValues, double delta)
        {
            Assert.AreEqual(expectedValues.Count, actualValues.Count);
            var dif = expectedValues.Zip(actualValues, (x, y) => Math.Abs(x - y)).Sum();
            var relDif = dif / expectedValues.Count;
            Assert.AreEqual(0, relDif, delta);
        }
    }
}