{
  description = "Nix flake for Tasks";
  inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs";

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        overrideAttrs.tasks = {
          src = lib.cleanSource ./.;
        };
      }
    );
}
