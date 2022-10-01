solution "ramanujan"
  architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

project "ramanujan"
  kind      "StaticLib"
  language  "C++"
  characterset ("MBCS")

  targetdir ("bin/" .. outputdir .. "/%{prj.name}")
  objdir ("bin-intermediate/" .. outputdir .. "/%{prj.name}")

  files
  {
    "source/**.h",
    "source/**.cpp"
  }

  includedirs
  {
    "source/includes"
  }

  filter { "files:**.c" }
    compileas "C++"

  filter "system:windows"
    cppdialect "C++17"
    staticruntime "On"
    systemversion "latest"

  filter "configurations:Debug"
	buildoptions "/MDd"
	symbols "On"

  filter "configurations:Release"
	buildoptions "/MD"
	optimize "On"
