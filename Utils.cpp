
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "Utils.h"


////////////////////////////////
//
// class Utils


/* quick sort function */

void Utils::quick_sort_min(int *ind, double *val, int size)
{
	int i, j;
	int t;
	double d;

	double pivot;

	if (size <= 24) {

		/* insertion sort */
		for (i=1;i<size;i++) {
			t = ind[i];
			d = val[i];
			for (j=i-1;j>=0;j--) {
				if (val[j] > d) {
					ind[j+1] = ind[j];
					val[j+1] = val[j];
				}
				else {
					break;
				}
			}
			ind[j+1] = t;
			val[j+1] = d;
		}
	}
	else {
		pivot = val[rand_int(size - 1)];

		i = -1;
		j = size;
		while (true) {
			do {
				i++;
			} while (val[i] < pivot);
			do {
				j--;
			} while (val[j] > pivot);
			if (i < j) {
				t = ind[i];
				d = val[i];
				ind[i] = ind[j];
				val[i] = val[j];
				ind[j] = t;
				val[j] = d;
			}
			else {
				if (i == 0) {
					i = 1;
				}
				break;
			}
		}

		quick_sort_min(ind, val, i);
		quick_sort_min(ind + i, val + i, size - i);
	}
}

void Utils::quick_sort_max(int *ind, double *val, int size)
{
	int i, j;
	int t;
	double d;

	double pivot;

	if (size <= 24) {

		/* insertion sort */
		for (i=1;i<size;i++) {
			t = ind[i];
			d = val[i];
			for (j=i-1;j>=0;j--) {
				if (val[j] < d) {
					ind[j+1] = ind[j];
					val[j+1] = val[j];
				}
				else {
					break;
				}
			}
			ind[j+1] = t;
			val[j+1] = d;
		}
	}
	else {
		pivot = val[rand_int(size - 1)];

		i = -1;
		j = size;
		while (true) {
			do {
				i++;
			} while (val[i] > pivot);
			do {
				j--;
			} while (val[j] < pivot);
			if (i < j) {
				t = ind[i];
				d = val[i];
				ind[i] = ind[j];
				val[i] = val[j];
				ind[j] = t;
				val[j] = d;
			}
			else {
				if (i == 0) {
					i = 1;
				}
				break;
			}
		}

		quick_sort_max(ind, val, i);
		quick_sort_max(ind + i, val + i, size - i);
	}
}

/* hash function */

#define hash_mix(a, b, c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
    b -= c; b -= a; b ^= (a<<8);  \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12); \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5);  \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}

unsigned long Utils::hash(register const unsigned char *key, register unsigned long len, register unsigned long init)
{
	register unsigned long a, b, c, l;
	/* set up the internal state */
	l = len;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = init;         /* the previous hash value */
	/* handle most of the key */
	while (l >= 12)
	{
		a += (key[0] + ((unsigned long) key[1]<<8) + ((unsigned long) key[2]<<16) + ((unsigned long) key[3]<<24));
		b += (key[4] + ((unsigned long) key[5]<<8) + ((unsigned long) key[6]<<16) + ((unsigned long) key[7]<<24));
		c += (key[8] + ((unsigned long) key[9]<<8) + ((unsigned long) key[10]<<16) + ((unsigned long) key[11]<<24));
		hash_mix(a, b, c);
		key += 12;
		l -= 12;
	}
	/* handle the last 11 bytes */
	c += len;
	switch(l)	/* all the case statements fall through */
	{
		case 11: c += ((unsigned long) key[10]<<24);
		case 10: c += ((unsigned long) key[9]<<16);
		case 9 : c += ((unsigned long) key[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b += ((unsigned long) key[7]<<24);
		case 7 : b += ((unsigned long) key[6]<<16);
		case 6 : b += ((unsigned long) key[5]<<8);
		case 5 : b += key[4];
		case 4 : a += ((unsigned long) key[3]<<24);
		case 3 : a += ((unsigned long) key[2]<<16);
		case 2 : a += ((unsigned long) key[1]<<8);
		case 1 : a += key[0];
		/* case 0: nothing left to add */
	}
	hash_mix(a,b,c);
	return c;
}


////////////////////////////////
//
// class Timer


int Timer::m_instance_num = 0;

void Timer::begin(const char *description)
{
	m_running = true;
	m_duration = 0;
	strcpy(m_description, description);
	Logger::debug("TIMER %d: %s ...", m_id, m_description);
	m_start = clock();
}

void Timer::end()
{
	pause();
	Logger::debug("TIMER %d: %s = %f seconds", m_id, m_description, m_duration);
}

void Timer::pause()
{
	if (m_running) {
		m_finish = clock();
		m_duration += (double) (m_finish - m_start) / CLOCKS_PER_SEC;
		m_running = false;
	}
}

void Timer::resume()
{
	if (!m_running) {
		m_start = clock();
		m_running = true;
	}
}


////////////////////////////////
//
// class Logger


bool Logger::m_logging = true;
int Logger::m_log_level = 3;
const int Logger::m_log_level_error = 1;
const int Logger::m_log_level_warning = 2;
const int Logger::m_log_level_info = 3;
const int Logger::m_log_level_verbose = 4;
const int Logger::m_log_level_debug = 5;
Timer Logger::m_timer[10];

void Logger::setLogLevel(int level)
{
	m_log_level = level;
}

void Logger::enableLogging()
{
	m_logging = true;
}

void Logger::disableLogging()
{
	m_logging = false;
}

void Logger::error(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_error, stderr, "ERROR: ", format, argptr);
	va_end(argptr);
}

void Logger::warning(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_warning, stderr, "WARNING: ", format, argptr);
	va_end(argptr);
}

void Logger::info(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_info, stdout, NULL, format, argptr);
	va_end(argptr);
}

void Logger::verbose(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_verbose, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::status(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_print(m_log_level_verbose, stderr, NULL, format, argptr);
	print(m_log_level_verbose, stderr, NULL, "               \r");
	va_end(argptr);
}

void Logger::debug(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_debug, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::print(int level, FILE *fp, const char *prompt, const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_print(level, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::println(int level, FILE *fp, const char *prompt, const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(level, stderr, NULL, format, argptr);
	va_end(argptr);
}

inline void Logger::_print(int level, FILE *fp, const char *prompt, const char *format, va_list argptr)
{
	if (m_log_level >= level) {
		if (prompt != NULL) fprintf(fp, prompt);
		vfprintf(fp, format, argptr);
	}
}

inline void Logger::_println(int level, FILE *fp, const char *prompt, const char *format, va_list argptr)
{
	if (m_log_level >= level) {
		if (prompt != NULL) fprintf(fp, prompt);
		vfprintf(fp, format, argptr);
		fprintf(fp, "\n");
	}
}

void Logger::beginTimer(int i, const char *description)
{
	m_timer[i].begin(description);
}

void Logger::endTimer(int i)
{
	m_timer[i].end();
}

void Logger::pauseTimer(int i)
{
	m_timer[i].pause();
}

void Logger::resumeTimer(int i)
{
	m_timer[i].resume();
}
