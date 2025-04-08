#pragma once

#include <cryptominisat.h>

using var_id = uint32_t;

inline CMSat::Lit pos( var_id id )
{
    return CMSat::Lit( id, false );
}

inline CMSat::Lit neg( var_id id )
{
	return CMSat::Lit( id, true );
}

