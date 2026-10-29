//! Native scalar configuration projection. Unsupported shapes fail compilation.
//! No business-type map or duplicated native field/default declarations.
#[macro_export]
macro_rules! project_parameters {
    (@resume {$name:ident} @adapter $(#[$doc:meta])* $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($(#[$doc])* $name, $native { $($(#[$attr])* $field: $ty),* });
    };
    (@resume {@flatten $(#[$doc:meta])* $name:ident, $native:path {$($parent:tt)*}, $field:ident} @$marker:ident $(#[$child_doc:meta])* $child:path { $($children:tt)* } methods {}) => {
        cosmolkit_tsify_macros::flatten_parameters!($(#[$doc])* $name,$native{$($parent)*},$field,{$($children)*});
    };
    (@resume {@flatten $(#[$doc:meta])* $name:ident, $native:path {$($parent:tt)*} factory, $field:ident} @$marker:ident $(#[$child_doc:meta])* $child:path { $($children:tt)* } methods {}) => {
        cosmolkit_tsify_macros::flatten_parameters!($(#[$doc])* $name,$native{$($parent)*} factory,$field,{$($children)*});
    };
    (@resume {$name:ident} @fields_factory $(#[$doc:meta])* $native:path { $($fields:tt)* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($(#[$doc])* $name,$native{$($fields)*} factory);
    };
    (@resume {$name:ident} @native_factory $(#[$doc:meta])* $native:path { $($fields:tt)* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($(#[$doc])* $name,$native{$($fields)*} factory);
    };
    (@resume {$name:ident} @native $(#[$doc:meta])* $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($(#[$doc])* $name, $native { $($(#[$attr])* $field: $ty),* });
    };
    (@resume {$name:ident} @fields $(#[$doc:meta])* $native:path { $($(#[$attr:meta])* $field:ident: $ty:ty),* } methods {}) => {
        cosmolkit_tsify_macros::parameters!($(#[$doc])* $name, $native { $($(#[$attr])* $field: $ty),* });
    };
}
