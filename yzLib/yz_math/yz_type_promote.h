/***********************************************************/
/**	\file
	\brief		type promote
	\details	type promotion using template according to 
				standard rules. 14 types in c++ standard

				bool < char < unsigened char < short < unsigned short < 
				int < unsigned int < long < unsigned long < long long < 
				unsigned long long < float < double < long double

				int is the minimal intermediate type, all types below int
				are promoted to int during calculation
	\author		Yizhong Zhang
	\date		5/20/2012
*/
/***********************************************************/
#ifndef __YZ_TYPE_PROMOTE_H__
#define __YZ_TYPE_PROMOTE_H__

namespace yz{

/**
	Type Promotion struct Macro: given two types, return the bigger one of them

	If one of the types is not one of 14 basic types, then return the type of T1

	type promote example, the return type is the bigger one of a and b:

	template<T1, T2>
	TYPE_PROMOTE(T1, T2) sum(T1 a, T2 b){
		return a + b;
	}
*/
#define TYPE_PROMOTE(T1, T2) typename yz::type_promote<T1, T2>::type

//	predefined macro to promote template T, T1, T2
//	these macros are partially used to overcome doxygen bugs
#define PROMOTE_T_TO_FLOAT		TYPE_PROMOTE(T, float)
#define PROMOTE_T1_T2			TYPE_PROMOTE(T1, T2)
#define PROMOTE_T1_T2_TO_FLOAT	TYPE_PROMOTE(PROMOTE_T1_T2, float)

/**
	Type Protion struct

	use TYPE_PROMOTE(T1, T2) for simplicity
*/
template<typename T1, typename T2>
struct type_promote{
	typedef T1 type;	///<	If T1 or T2 is not basic type, we don't know how to promote, so just use type T1
};

/**
	set type promotion rules

	in this file, we set promotion rules for basic data types, 
	if we want to promote other types of data, call
	SET_PROMOTION_RULES(T1, T2, promoted_type) after the types
	are defined. Then the rules will be set
*/
#define SET_PROMOTION_RULES(T1, T2, promoted_type) \
template<>struct type_promote<T1, T2>	{typedef promoted_type type;};

SET_PROMOTION_RULES(bool,				bool,				int					)
SET_PROMOTION_RULES(bool,				char,				int					)
SET_PROMOTION_RULES(bool,				unsigned char,		int					)
SET_PROMOTION_RULES(bool,				short,				int					)
SET_PROMOTION_RULES(bool,				unsigned short,		int					)
SET_PROMOTION_RULES(bool,				int,				int					)
SET_PROMOTION_RULES(bool,				unsigned int,		unsigned int		)
SET_PROMOTION_RULES(bool,				long,				long				)
SET_PROMOTION_RULES(bool,				unsigned long,		unsigned long		)
SET_PROMOTION_RULES(bool,				long long,			long long			)
SET_PROMOTION_RULES(bool,				unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(bool,				float,				float				)
SET_PROMOTION_RULES(bool,				double,				double				)
SET_PROMOTION_RULES(bool,				long double,		long double			)

SET_PROMOTION_RULES(char,				bool,				int					)
SET_PROMOTION_RULES(char,				char,				int					)
SET_PROMOTION_RULES(char,				unsigned char,		int					)
SET_PROMOTION_RULES(char,				short,				int					)
SET_PROMOTION_RULES(char,				unsigned short,		int					)
SET_PROMOTION_RULES(char,				int,				int					)
SET_PROMOTION_RULES(char,				unsigned int,		unsigned int		)
SET_PROMOTION_RULES(char,				long,				long				)
SET_PROMOTION_RULES(char,				unsigned long,		unsigned long		)
SET_PROMOTION_RULES(char,				long long,			long long			)
SET_PROMOTION_RULES(char,				unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(char,				float,				float				)
SET_PROMOTION_RULES(char,				double,				double				)
SET_PROMOTION_RULES(char,				long double,		long double			)

SET_PROMOTION_RULES(unsigned char,		bool,				int					)
SET_PROMOTION_RULES(unsigned char,		char,				int					)
SET_PROMOTION_RULES(unsigned char,		unsigned char,		int					)
SET_PROMOTION_RULES(unsigned char,		short,				int					)
SET_PROMOTION_RULES(unsigned char,		unsigned short,		int					)
SET_PROMOTION_RULES(unsigned char,		int,				int					)
SET_PROMOTION_RULES(unsigned char,		unsigned int,		unsigned int		)
SET_PROMOTION_RULES(unsigned char,		long,				long				)
SET_PROMOTION_RULES(unsigned char,		unsigned long,		unsigned long		)
SET_PROMOTION_RULES(unsigned char,		long long,			long long			)
SET_PROMOTION_RULES(unsigned char,		unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(unsigned char,		float,				float				)
SET_PROMOTION_RULES(unsigned char,		double,				double				)
SET_PROMOTION_RULES(unsigned char,		long double,		long double			)

SET_PROMOTION_RULES(short,				bool,				int					)
SET_PROMOTION_RULES(short,				char,				int					)
SET_PROMOTION_RULES(short,				unsigned char,		int					)
SET_PROMOTION_RULES(short,				short,				int					)
SET_PROMOTION_RULES(short,				unsigned short,		int					)
SET_PROMOTION_RULES(short,				int,				int					)
SET_PROMOTION_RULES(short,				unsigned int,		unsigned int		)
SET_PROMOTION_RULES(short,				long,				long				)
SET_PROMOTION_RULES(short,				unsigned long,		unsigned long		)
SET_PROMOTION_RULES(short,				long long,			long long			)
SET_PROMOTION_RULES(short,				unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(short,				float,				float				)
SET_PROMOTION_RULES(short,				double,				double				)
SET_PROMOTION_RULES(short,				long double,		long double			)

SET_PROMOTION_RULES(unsigned short,		bool,				int					)
SET_PROMOTION_RULES(unsigned short,		char,				int					)
SET_PROMOTION_RULES(unsigned short,		unsigned char,		int					)
SET_PROMOTION_RULES(unsigned short,		short,				int					)
SET_PROMOTION_RULES(unsigned short,		unsigned short,		int					)
SET_PROMOTION_RULES(unsigned short,		int,				int					)
SET_PROMOTION_RULES(unsigned short,		unsigned int,		unsigned int		)
SET_PROMOTION_RULES(unsigned short,		long,				long				)
SET_PROMOTION_RULES(unsigned short,		unsigned long,		unsigned long		)
SET_PROMOTION_RULES(unsigned short,		long long,			long long			)
SET_PROMOTION_RULES(unsigned short,		unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(unsigned short,		float,				float				)
SET_PROMOTION_RULES(unsigned short,		double,				double				)
SET_PROMOTION_RULES(unsigned short,		long double,		long double			)

SET_PROMOTION_RULES(int,				bool,				int					)
SET_PROMOTION_RULES(int,				char,				int					)
SET_PROMOTION_RULES(int,				unsigned char,		int					)
SET_PROMOTION_RULES(int,				short,				int					)
SET_PROMOTION_RULES(int,				unsigned short,		int					)
SET_PROMOTION_RULES(int,				int,				int					)
SET_PROMOTION_RULES(int,				unsigned int,		unsigned int		)
SET_PROMOTION_RULES(int,				long,				long				)
SET_PROMOTION_RULES(int,				unsigned long,		unsigned long		)
SET_PROMOTION_RULES(int,				long long,			long long			)
SET_PROMOTION_RULES(int,				unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(int,				float,				float				)
SET_PROMOTION_RULES(int,				double,				double				)
SET_PROMOTION_RULES(int,				long double,		long double			)

SET_PROMOTION_RULES(unsigned int,		bool,				unsigned int		)
SET_PROMOTION_RULES(unsigned int,		char,				unsigned int		)
SET_PROMOTION_RULES(unsigned int,		unsigned char,		unsigned int		)
SET_PROMOTION_RULES(unsigned int,		short,				unsigned int		)
SET_PROMOTION_RULES(unsigned int,		unsigned short,		unsigned int		)
SET_PROMOTION_RULES(unsigned int,		int,				unsigned int		)
SET_PROMOTION_RULES(unsigned int,		unsigned int,		unsigned int		)
SET_PROMOTION_RULES(unsigned int,		long,				long				)
SET_PROMOTION_RULES(unsigned int,		unsigned long,		unsigned long		)
SET_PROMOTION_RULES(unsigned int,		long long,			long long			)
SET_PROMOTION_RULES(unsigned int,		unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(unsigned int,		float,				float				)
SET_PROMOTION_RULES(unsigned int,		double,				double				)
SET_PROMOTION_RULES(unsigned int,		long double,		long double			)

SET_PROMOTION_RULES(long,				bool,				long				)
SET_PROMOTION_RULES(long,				char,				long				)
SET_PROMOTION_RULES(long,				unsigned char,		long				)
SET_PROMOTION_RULES(long,				short,				long				)
SET_PROMOTION_RULES(long,				unsigned short,		long				)
SET_PROMOTION_RULES(long,				int,				long				)
SET_PROMOTION_RULES(long,				unsigned int,		long				)
SET_PROMOTION_RULES(long,				long,				long				)
SET_PROMOTION_RULES(long,				unsigned long,		unsigned long		)
SET_PROMOTION_RULES(long,				long long,			long long			)
SET_PROMOTION_RULES(long,				unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(long,				float,				float				)
SET_PROMOTION_RULES(long,				double,				double				)
SET_PROMOTION_RULES(long,				long double,		long double			)

SET_PROMOTION_RULES(unsigned long,		bool,				unsigned long		)
SET_PROMOTION_RULES(unsigned long,		char,				unsigned long		)
SET_PROMOTION_RULES(unsigned long,		unsigned char,		unsigned long		)
SET_PROMOTION_RULES(unsigned long,		short,				unsigned long		)
SET_PROMOTION_RULES(unsigned long,		unsigned short,		unsigned long		)
SET_PROMOTION_RULES(unsigned long,		int,				unsigned long		)
SET_PROMOTION_RULES(unsigned long,		unsigned int,		unsigned long		)
SET_PROMOTION_RULES(unsigned long,		long,				unsigned long		)
SET_PROMOTION_RULES(unsigned long,		unsigned long,		unsigned long		)
SET_PROMOTION_RULES(unsigned long,		long long,			long long			)
SET_PROMOTION_RULES(unsigned long,		unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(unsigned long,		float,				float				)
SET_PROMOTION_RULES(unsigned long,		double,				double				)
SET_PROMOTION_RULES(unsigned long,		long double,		long double			)

SET_PROMOTION_RULES(long long,			bool,				long long			)
SET_PROMOTION_RULES(long long,			char,				long long			)
SET_PROMOTION_RULES(long long,			unsigned char,		long long			)
SET_PROMOTION_RULES(long long,			short,				long long			)
SET_PROMOTION_RULES(long long,			unsigned short,		long long			)
SET_PROMOTION_RULES(long long,			int,				long long			)
SET_PROMOTION_RULES(long long,			unsigned int,		long long			)
SET_PROMOTION_RULES(long long,			long,				long long			)
SET_PROMOTION_RULES(long long,			unsigned long,		long long			)
SET_PROMOTION_RULES(long long,			long long,			long long			)
SET_PROMOTION_RULES(long long,			unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(long long,			float,				float				)
SET_PROMOTION_RULES(long long,			double,				double				)
SET_PROMOTION_RULES(long long,			long double,		long double			)

SET_PROMOTION_RULES(unsigned long long,	bool,				unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	char,				unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	unsigned char,		unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	short,				unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	unsigned short,		unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	int,				unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	unsigned int,		unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	long,				unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	unsigned long,		unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	long long,			unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	unsigned long long,	unsigned long long	)
SET_PROMOTION_RULES(unsigned long long,	float,				float				)
SET_PROMOTION_RULES(unsigned long long,	double,				double				)
SET_PROMOTION_RULES(unsigned long long,	long double,		long double			)

SET_PROMOTION_RULES(float,				bool,				float				)
SET_PROMOTION_RULES(float,				char,				float				)
SET_PROMOTION_RULES(float,				unsigned char,		float				)
SET_PROMOTION_RULES(float,				short,				float				)
SET_PROMOTION_RULES(float,				unsigned short,		float				)
SET_PROMOTION_RULES(float,				int,				float				)
SET_PROMOTION_RULES(float,				unsigned int,		float				)
SET_PROMOTION_RULES(float,				long,				float				)
SET_PROMOTION_RULES(float,				unsigned long,		float				)
SET_PROMOTION_RULES(float,				long long,			float				)
SET_PROMOTION_RULES(float,				unsigned long long,	float				)
SET_PROMOTION_RULES(float,				float,				float				)
SET_PROMOTION_RULES(float,				double,				double				)
SET_PROMOTION_RULES(float,				long double,		long double			)

SET_PROMOTION_RULES(double,				bool,				double				)
SET_PROMOTION_RULES(double,				char,				double				)
SET_PROMOTION_RULES(double,				unsigned char,		double				)
SET_PROMOTION_RULES(double,				short,				double				)
SET_PROMOTION_RULES(double,				unsigned short,		double				)
SET_PROMOTION_RULES(double,				int,				double				)
SET_PROMOTION_RULES(double,				unsigned int,		double				)
SET_PROMOTION_RULES(double,				long,				double				)
SET_PROMOTION_RULES(double,				unsigned long,		double				)
SET_PROMOTION_RULES(double,				long long,			double				)
SET_PROMOTION_RULES(double,				unsigned long long,	double				)
SET_PROMOTION_RULES(double,				float,				double				)
SET_PROMOTION_RULES(double,				double,				double				)
SET_PROMOTION_RULES(double,				long double,		long double			)

SET_PROMOTION_RULES(long double,		bool,				long double			)
SET_PROMOTION_RULES(long double,		char,				long double			)
SET_PROMOTION_RULES(long double,		unsigned char,		long double			)
SET_PROMOTION_RULES(long double,		short,				long double			)
SET_PROMOTION_RULES(long double,		unsigned short,		long double			)
SET_PROMOTION_RULES(long double,		int,				long double			)
SET_PROMOTION_RULES(long double,		unsigned int,		long double			)
SET_PROMOTION_RULES(long double,		long,				long double			)
SET_PROMOTION_RULES(long double,		unsigned long,		long double			)
SET_PROMOTION_RULES(long double,		long long,			long double			)
SET_PROMOTION_RULES(long double,		unsigned long long,	long double			)
SET_PROMOTION_RULES(long double,		float,				long double			)
SET_PROMOTION_RULES(long double,		double,				long double			)
SET_PROMOTION_RULES(long double,		long double,		long double			)

#undef SET_PROMOTION_RULES

/**
	Check Whether Type Is bool

	to use: Is_bool<T>::check_type
*/
template<typename T>
struct Is_bool{
	static const bool check_type = false;
};
template<> struct Is_bool<bool>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is char

	to use: Is_char<T>::check_type
*/
template<typename T>
struct Is_char{
	static const bool check_type = false;
};
template<> struct Is_char<char>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is unsigned char

	to use: Is_unsigned_char<T>::check_type
*/
template<typename T>
struct Is_unsigned_char{
	static const bool check_type = false;
};
template<> struct Is_unsigned_char<unsigned char>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is short

	to use: Is_short<T>::check_type
*/
template<typename T>
struct Is_short{
	static const bool check_type = false;
};
template<> struct Is_short<short>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is unsigned_short

	to use: Is_unsigned_short<T>::check_type
*/
template<typename T>
struct Is_unsigned_short{
	static const bool check_type = false;
};
template<> struct Is_unsigned_short<unsigned short>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is int

	to use: Is_int<T>::check_type
*/
template<typename T>
struct Is_int{
	static const bool check_type = false;
};
template<> struct Is_int<int>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is unsigned int

	to use: Is_unsigned_int<T>::check_type
*/
template<typename T>
struct Is_unsigned_int{
	static const bool check_type = false;
};
template<> struct Is_unsigned_int<unsigned int>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is long

	to use: Is_long<T>::check_type
*/
template<typename T>
struct Is_long{
	static const bool check_type = false;
};
template<> struct Is_long<long>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is unsigned long

	to use: Is_unsigned_long<T>::check_type
*/
template<typename T>
struct Is_unsigned_long{
	static const bool check_type = false;
};
template<> struct Is_unsigned_long<unsigned long>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is long long

	to use: Is_long_long<T>::check_type
*/
template<typename T>
struct Is_long_long{
	static const bool check_type = false;
};
template<> struct Is_long_long<long long>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is unsigned long long

	to use: Is_unsigned_long_long<T>::check_type
*/
template<typename T>
struct Is_unsigned_long_long{
	static const bool check_type = false;
};
template<> struct Is_unsigned_long_long<unsigned long long>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is float

	to use: Is_float<T>::check_type
*/
template<typename T>
struct Is_float{
	static const bool check_type = false;
};
template<> struct Is_float<float>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is double

	to use: Is_double<T>::check_type
*/
template<typename T>
struct Is_double{
	static const bool check_type = false;
};
template<> struct Is_double<double>{
	static const bool check_type = true;
};

/**
	Check Whether Type Is long double

	to use: Is_long_double<T>::check_type
*/
template<typename T>
struct Is_long_double{
	static const bool check_type = false;
};
template<> struct Is_long_double<long double>{
	static const bool check_type = true;
};


}	//	namespace yz

#endif	//	__YZ_TYPE_PROMOTE_H__