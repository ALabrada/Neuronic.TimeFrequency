namespace Neuronic.TimeFrequency
{
    public interface IStopCriteria<T>
    {
        bool ShouldStop(T state);
    }
}