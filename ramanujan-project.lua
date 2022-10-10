project "ramanujan"
kind "StaticLib"
language "C++"
characterset("MBCS")

targetdir("bin/" .. outputdir .. "/%{prj.name}")
objdir("bin-intermediate/" .. outputdir .. "/%{prj.name}")

files
{
  "source/**.h",
  "source/**.cpp"
}

includedirs
{
  "source"
}

filter { "files:**.c" }
compileas "C++"

filter "system:windows"
cppdialect "C++17"
-- staticruntime "On"
systemversion "latest"

filter "configurations:Debug"
symbols "On"
staticruntime "off"
runtime "Debug"

filter "configurations:Release"
optimize "On"
staticruntime "off"
runtime "Release"