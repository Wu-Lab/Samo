
#ifndef __UTILS_H
#define __UTILS_H


#include <ctime>
#include <vector>
#include <iosfwd>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <boost/shared_ptr.hpp>

#include "Matrix.h"


using namespace std;


namespace std {
	namespace tr1 = ::boost;
}


class Utils {
public:
	template <class T> static int rand_int(T max) { return (int) ((double) rand() * max / RAND_MAX); }

	template <class T> static T min (T i, T j) { return ((i < j) ? i : j); }
	template <class T> static T max (T i, T j) { return ((i > j) ? i : j); }

	template <class T> static void swap(T &i, T &j) { T temp = i; i = j; j = temp; }

	static void quick_sort_min(int *ind, double *val, int size);
	static void quick_sort_max(int *ind, double *val, int size);

	static unsigned long hash(const char *key) { return hash((const unsigned char *) key, strlen(key)); }
	static unsigned long hash(register const unsigned char *key, register unsigned long len, register unsigned long init = 0);
	static unsigned long hash_size(int n) { return (unsigned long) 1 << n; }
	static unsigned long hash_mask(int n) { return hash_size(n)-1; }
};


class Timer {
	static int m_instance_num;
	int m_id;
	clock_t m_start, m_finish;
	bool m_running;
	double m_duration;
	char m_description[80];

public:
	Timer() { m_id = m_instance_num++; }
	void begin(const char *description);
	void end();
	void pause();
	void resume();
	double time() { return m_duration; }
};


class Logger {
	static bool m_logging;
	static int m_log_level;
	static const int m_log_level_error;
	static const int m_log_level_warning;
	static const int m_log_level_info;
	static const int m_log_level_verbose;
	static const int m_log_level_debug;

	static Timer m_timer[10];

public:
	static int log_level() { return m_log_level; }
	static int log_level_error() { return m_log_level_error; }
	static int log_level_warning() { return m_log_level_warning; }
	static int log_level_info() { return m_log_level_info; }
	static int log_level_verbose() { return m_log_level_verbose; }
	static int log_level_debug() { return m_log_level_debug; }

	static Timer &timer(int i) { return m_timer[i]; }

	static void error(const char *format, ...);
	static void warning(const char *format, ...);
	static void info(const char *format, ...);
	static void verbose(const char *format, ...);
	static void status(const char *format, ...);
	static void debug(const char *format, ...);

	static void print(int level, FILE *fp, const char *prompt, const char *format, ...);
	static void println(int level, FILE *fp, const char *prompt, const char *format, ...);

	static void setLogLevel(int level);
	static void enableLogging();
	static void disableLogging();

	static bool isDebug() { return (m_log_level >= m_log_level_debug); }

	static void beginTimer(int i, const char *description);
	static void endTimer(int i);
	static void pauseTimer(int i);
	static void resumeTimer(int i);

private:
	static void _print(int level, FILE *fp, const char *prompt, const char *format, va_list argptr);
	static void _println(int level, FILE *fp, const char *prompt, const char *format, va_list argptr);
};


template<typename T>
ostream &operator<<(ostream &os, const vector<T> &v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " ")); 
    return os;
};


// functions for string

void string_replace(string &str, const string &src, const string &dst);
string int2str(int num);
int str2int(const string &str);


#endif // __UTILS_H
