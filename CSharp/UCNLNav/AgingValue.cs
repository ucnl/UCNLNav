using System;
using System.Text;

namespace UCNLNav
{
    public interface IAging
    {
        bool IsInitialized { get; }
        bool IsObsolete { get; }
    }

    
    public class AgingValue<T> : IAging
    {
        #region Properties

        public static readonly int DefaultMinSecToShowAge = 5;
        public static readonly int DefaultObsoleteIntervalSec = 300;

        T val;
        public T Value
        {
            get { return val; }
            set
            {
                val = value;
                IsInitialized = true;
                TimeStamp = DateTime.Now;
            }
        }

        public string AccessTag { get; set; }

        public bool IsInitialized { get; private set; }

        public DateTime TimeStamp { get; private set; }

        public int ObsoleteIntervalSec { get; set; }

        public void ForceUpdate()
        {
            TimeStamp = DateTime.Now;
        }

        public int MinSecToShowAge { get; set; }

        public TimeSpan Age
        {
            get
            {
                return DateTime.Now.Subtract(TimeStamp);
            }
        }

        public bool IsObsolete
        {
            get
            {
                return Age.TotalSeconds > ObsoleteIntervalSec;
            }
        }

        public bool IsInitializedAndNotObsolete
        {
            get { return IsInitialized && !IsObsolete; }
        }

        Func<T, string> CustomFormatter;        

        #endregion

        #region Constructor

        public AgingValue(int minTimeToShowAgeSec, int obsoleteIntervalSec, Func<T, string> customFormatter)
        {
            if (customFormatter == null)
                throw new ArgumentNullException("customFormatter");

            IsInitialized = false;
            MinSecToShowAge = minTimeToShowAgeSec;
            ObsoleteIntervalSec = obsoleteIntervalSec;
            CustomFormatter = customFormatter;
        }

        public AgingValue()
            : this(DefaultMinSecToShowAge, DefaultObsoleteIntervalSec, ((v) => { return v.ToString(); }))
        {
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            if (!IsInitialized)
            {
                return "- - -";
            }
            else if (!string.IsNullOrEmpty(AccessTag))
            {
                return AccessTag;
            }
            else
            {
                StringBuilder sb = new StringBuilder();
                sb.Append(CustomFormatter(val));

                TimeSpan age = Age;

                if (Age.TotalSeconds > MinSecToShowAge)
                {
                    if (age.TotalSeconds > ObsoleteIntervalSec)
                        sb.Append(" (OBS)");
                    else
                        sb.AppendFormat(" ({0:00}:{1:00})", age.Minutes, age.Seconds);
                }

                return sb.ToString();
            }
        }

        #endregion
    }
}
