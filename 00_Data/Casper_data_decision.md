# Casper dataset decision needed

The supplied discharge file contains eight South Fork Caspar Creek watersheds:
`OGI`, `POR`, `RIC`, `SEQ`, `TRE`, `UQL`, `WIL`, and `ZIE`.

The supplied precipitation file contains two gauges: `S620` and `S640`.
Neither the local files nor the project notes provide a watershed-to-gauge mapping.
The TE workflow therefore does not run Casper until the precipitation rule is chosen
(one gauge for every watershed, the mean of both gauges, or an explicit mapping).

