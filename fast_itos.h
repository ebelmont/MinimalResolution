#pragma once
#include<string>
#include<cstdlib>

//fast integer-to-string conversion avoiding std::to_string's vsnprintf-based
//implementation (locale-aware formatted printing, much slower than plain
//digit extraction) in this libstdc++ build.
inline std::string fast_itos(long long v){
	if(v==0) return std::string("0");
	bool neg = v<0;
	unsigned long long uv = neg ? (unsigned long long)(-(v+1))+1ULL : (unsigned long long)v;
	char buf[24];
	int pos = 24;
	while(uv){
		buf[--pos] = '0' + (char)(uv%10);
		uv/=10;
	}
	if(neg) buf[--pos] = '-';
	return std::string(buf+pos, 24-pos);
}

inline std::string fast_itos(int v){ return fast_itos((long long)v); }
inline std::string fast_itos(unsigned v){ return fast_itos((long long)v); }
