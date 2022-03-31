/***********************************************************/
/**	\file
	\brief		HASH
	\author		Yizhong Zhang
	\date		11/16/2017
*/
/***********************************************************/
#ifndef __YZ_HASH_H__
#define __YZ_HASH_H__

#include <stdint.h>

namespace yz{ namespace utils{

//	========================================
///@{
/**	@name HASH related
*/
//	========================================

/**
	Bitwise hasher for stl, can be used for arbitrary data type

	For example: \n
	std::unordered_map<yz::Vec3d, int, BitwiseHasher<yz::Vec3d>> vec3d_hash;
*/
template <class T>
class BitwiseHasher
{
public:
	/**
	From vs2017, std::_Hash_seq() is no longer available. so use this function instead
	*/
	inline size_t fnv1a_hash_bytes(const unsigned char * first, size_t count) const{
#if UINTPTR_MAX == 0xFFFFFFFFFFFFFFFF	//	64 bit
		static_assert(sizeof(size_t) == 8, "This code is for 64-bit size_t.");
		const size_t fnv_offset_basis = 14695981039346656037ULL;
		const size_t fnv_prime = 1099511628211ULL;

#elif UINTPTR_MAX == 0xFFFFFFFF			//	32 bit
		static_assert(sizeof(size_t) == 4, "This code is for 32-bit size_t.");
		const size_t fnv_offset_basis = 2166136261U;
		const size_t fnv_prime = 16777619U;

#else
		#error yz::utils::BitwiseHasher is only defined in 64 and 32 bits

#endif 

		size_t result = fnv_offset_basis;
		for (size_t next = 0; next < count; ++next)
		{	// fold in another byte
			result ^= (size_t)first[next];
			result *= fnv_prime;
		}
		return (result);
	}

	std::size_t operator()(const T& k) const
	{
		return fnv1a_hash_bytes((const unsigned char *)&k, sizeof(T));
	}
};


///@}

}}	//	end namespace yz::utils

#endif	//	__YZ_HASH_H__